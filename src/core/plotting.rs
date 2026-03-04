// src/solver/plotting.rs

use plotters::prelude::*;
use std::error::Error;
use std::path::Path;

pub fn plot_mr_curve(
    radii_km: &[f64],
    masses_msun: &[f64],
    output_path: &str,
) -> Result<(), Box<dyn Error>> {
    
    let path = Path::new(output_path);
    if let Some(parent) = path.parent() {
       std::fs::create_dir_all(parent)?;
    }

    // --- High DPI Configuration ---
    // Scale factor of 3 turns an 800x600 plot into a crisp 2400x1800 plot
    let scale = 1; 
    let width = 800 * scale;
    let height = 600 * scale;

    let root = SVGBackend::new(output_path, (width, height)).into_drawing_area();
    // let root = BitMapBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let r_min = 8.0; 
    let r_max = 14.0;
     // let r_max = radii_km.iter().fold(0./0., |a: f64, b| a.max(*b)) * 1.05;
    let m_max = masses_msun.iter().fold(0./0., |a: f64, b| a.max(*b)) * 1.05;

    // Multiply all font sizes and margins by the 'scale' variable
    let mut chart = ChartBuilder::on(&root)
        .caption("Mass-Radius Relation", ("sans-serif", 24 * scale).into_font())
        .margin(15 * scale)
        .x_label_area_size(40 * scale)
        .y_label_area_size(50 * scale)
        .build_cartesian_2d(r_min..r_max, 0.0f64..m_max)?;

    chart
        .configure_mesh()
        .x_desc("Radius [km]")
        // \u{2299} is the Unicode character for the Sun symbol ⊙
        // \u{03C1} would be the Greek letter rho ρ, etc.
        .y_desc("Mass [M\u{2299}]") 
        .axis_desc_style(("sans-serif", 16 * scale))
        .draw()?;

    let plot_blue = RGBColor(31, 119, 180);

    let data_iter = radii_km.iter().zip(masses_msun.iter());

    chart
        .draw_series(LineSeries::new(
            data_iter.map(|(r, m)| (*r, *m)),
            plot_blue.stroke_width(2 * scale), // Scale the line thickness
        ))?
        .label("EoS Model")
        .legend(move |(x, y)| PathElement::new(
            vec![(x, y), (x + (20 * scale) as i32, y)], 
            plot_blue.stroke_width(2 * scale)
        ));

    chart
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .label_font(("sans-serif", 12 * scale)) // Scale the legend font
        .draw()?;

    root.present()?;
    Ok(())
}


#[derive(Clone)]
pub struct CurveData {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub label: String,
}

pub struct Artist {
    output_path: String,
    title: String,
    x_label: String,
    y_label: String,
    x_min: Option<f64>,
    x_max: Option<f64>,
    y_min: Option<f64>,
    y_max: Option<f64>,
    curves: Vec<CurveData>,
}

impl Artist {
    pub fn new(output_path: &str, title: &str) -> Self {
        Self {
            output_path: output_path.to_string(),
            title: title.to_string(),
            x_label: "Radius [km]".to_string(),
            y_label: "Mass [M\u{2299}]".to_string(),
            x_min: None,
            x_max: None,
            y_min: None,
            y_max: None,
            curves: Vec::new(),
        }
    }

    pub fn with_x_label(mut self, label: &str) -> Self {
        self.x_label = label.to_string();
        self
    }

    pub fn with_y_label(mut self, label: &str) -> Self {
        self.y_label = label.to_string();
        self
    }

    /// Ativa o ajuste automático (Limpa limites manuais)
    pub fn autoscale(mut self) -> Self {
        self.x_min = None;
        self.x_max = None;
        self.y_min = None;
        self.y_max = None;
        self
    }

    pub fn with_x_range(mut self, min: f64, max: f64) -> Self {
        self.x_min = Some(min);
        self.x_max = Some(max);
        self
    }

    pub fn add_curve(mut self, x: &[f64], y: &[f64], label: &str) -> Self {
        self.curves.push(CurveData {
            x: x.to_vec(),
            y: y.to_vec(),
            label: label.to_string(),
        });
        self
    }

    pub fn plot(&self) -> Result<(), Box<dyn Error>> {
        if self.curves.is_empty() { return Err("Sem curvas.".into()); }

        let path = Path::new(&self.output_path);
        if let Some(parent) = path.parent() { std::fs::create_dir_all(parent)?; }

        // --- LÓGICA DE ESCALA INTELIGENTE ---
        let mut d_x_min = f64::MAX;
        let mut d_x_max = f64::MIN;
        let mut d_y_min = f64::MAX;
        let mut d_y_max = f64::MIN;

        for curve in &self.curves {
            for (&xi, &yi) in curve.x.iter().zip(curve.y.iter()) {
                if xi.is_finite() && yi.is_finite() {
                    // FILTRO DE OUTLIERS: Ignora raios > 50km no MR para não distorcer o plot
                    if self.x_label.contains("Radius") && xi > 50.0 { continue; }
                    
                    if xi < d_x_min { d_x_min = xi; }
                    if xi > d_x_max { d_x_max = xi; }
                    if yi < d_y_min { d_y_min = yi; }
                    if yi > d_y_max { d_y_max = yi; }
                }
            }
        }

        // Se os dados estiverem vazios ou inválidos
        if d_x_min == f64::MAX { d_x_min = 0.0; d_x_max = 1.0; d_y_max = 1.0; }

        // Definição final dos limites (Prioridade para manual, senão automático com folga)
        let x_start = self.x_min.unwrap_or(if self.x_label.contains("Radius") { d_x_min.max(8.0) } else { d_x_min });
        let x_end = self.x_max.unwrap_or(d_x_max * 1.02);
        let y_start = self.y_min.unwrap_or(0.0); // Massa e Pressão geralmente começam em 0
        let y_end = self.y_max.unwrap_or(d_y_max * 1.05);

        // --- RENDERIZAÇÃO ---
        let root = SVGBackend::new(&self.output_path, (800, 600)).into_drawing_area();
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .caption(&self.title, ("sans-serif", 20))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(50)
            .build_cartesian_2d(x_start..x_end, y_start..y_end)?;

        chart.configure_mesh()
            .x_desc(&self.x_label).y_desc(&self.y_label)
            .light_line_style(&WHITE.mix(0.1)) // Deixa o grid mais suave
            .draw()?;

        let colors = [&BLUE, &RED, &GREEN, &MAGENTA, &CYAN, &BLACK];
        for (i, curve) in self.curves.iter().enumerate() {
            let color = colors[i % colors.len()];
            // Downsampling para performance
            let step = (curve.x.len() / 2000).max(1);
            let data = curve.x.iter().zip(curve.y.iter()).enumerate()
                .filter(|(idx, (x, y))| idx % step == 0 && x.is_finite() && y.is_finite())
                .map(|(_, (&x, &y))| (x, y));

            chart.draw_series(LineSeries::new(data, color.stroke_width(2)))?
                .label(&curve.label)
                .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2)));
        }

        chart.configure_series_labels().position(SeriesLabelPosition::UpperRight).draw()?;
        root.present()?;
        Ok(())
    }
}