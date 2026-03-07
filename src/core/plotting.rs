// src/solver/plotting.rs

use plotters::prelude::*;
use std::error::Error;
use std::path::Path;

// ============================================================================
// FUNÇÕES INDIVIDUAIS DE PLOTAGEM (Única Curva)
// ============================================================================

/// Plota uma curva Massa-Raio padrão (Escala Linear)
pub fn plot_mr_curve(
    radii_km: &[f64],
    masses_msun: &[f64],
    output_path: &str,
) -> Result<(), Box<dyn Error>> {
    let path = Path::new(output_path);
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let scale = 1; 
    let width = 800 * scale;
    let height = 600 * scale;

    let root = SVGBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let r_min = 8.0; 
    let r_max = 14.0;
    let m_max = masses_msun.iter().fold(0./0., |a: f64, b| a.max(*b)) * 1.05;

    let mut chart = ChartBuilder::on(&root)
        .caption("Mass-Radius Relation", ("sans-serif", 24 * scale).into_font())
        .margin(15 * scale)
        .x_label_area_size(40 * scale)
        .y_label_area_size(50 * scale)
        .build_cartesian_2d(r_min..r_max, 0.0f64..m_max)?;

    chart.configure_mesh()
        .x_desc("Radius [km]")
        .y_desc("Mass [M\u{2299}]") 
        .axis_desc_style(("sans-serif", 16 * scale))
        .draw()?;

    let plot_blue = RGBColor(31, 119, 180);
    let data_iter = radii_km.iter().zip(masses_msun.iter());

    chart.draw_series(LineSeries::new(
        data_iter.map(|(r, m)| (*r, *m)),
        plot_blue.stroke_width(2 * scale),
    ))?
    .label("MR Model")
    .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + (20 * scale) as i32, y)], plot_blue.stroke_width(2 * scale)));

    chart.configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .label_font(("sans-serif", 12 * scale))
        .draw()?;

    root.present()?;
    Ok(())
}

/// Plota uma curva de Equação de Estado (EoS) com suporte opcional a Log-Log
pub fn plot_eos_curve(
    eps_array: &[f64],
    p_array: &[f64],
    output_path: &str,
    use_logscale: bool,
) -> Result<(), Box<dyn Error>> {
    let path = Path::new(output_path);
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let scale = 1;
    let width = 800 * scale;
    let height = 600 * scale;
    let root = SVGBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut d_x_min = f64::MAX; let mut d_x_max = f64::MIN;
    let mut d_y_min = f64::MAX; let mut d_y_max = f64::MIN;

    let mut clean_data = Vec::new();

    // Filtro de dados
    for (&e, &p) in eps_array.iter().zip(p_array.iter()) {
        if e.is_finite() && p.is_finite() {
            if use_logscale && (e <= 0.0 || p <= 0.0) { continue; } // Log não aceita <= 0
            
            clean_data.push((e, p));
            if e < d_x_min { d_x_min = e; }
            if e > d_x_max { d_x_max = e; }
            if p < d_y_min { d_y_min = p; }
            if p > d_y_max { d_y_max = p; }
        }
    }

    if clean_data.is_empty() { return Err("Nenhum dado válido para plotar.".into()); }

    let plot_red = RGBColor(214, 39, 40);

    // Renderização separada devido aos tipos genéricos do plotters para Log Coord
    if use_logscale {
        let mut chart = ChartBuilder::on(&root)
            .caption("Equation of State", ("sans-serif", 24 * scale))
            .margin(15 * scale)
            .x_label_area_size(40 * scale)
            .y_label_area_size(50 * scale)
            .build_cartesian_2d((d_x_min * 0.8..d_x_max * 1.2).log_scale(), (d_y_min * 0.8..d_y_max * 1.5).log_scale())?;

        chart.configure_mesh()
            .x_desc("Energy Density \u{3B5} [MeV/fm\u{00B3}]")
            .y_desc("Pressure P [MeV/fm\u{00B3}]")
            .draw()?;

        chart.draw_series(LineSeries::new(clean_data, plot_red.stroke_width(2 * scale)))?
            .label("EoS Model")
            .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], plot_red.stroke_width(2)));

        chart.configure_series_labels().position(SeriesLabelPosition::UpperLeft).draw()?;
    } else {
        let mut chart = ChartBuilder::on(&root)
            .caption("Equation of State", ("sans-serif", 24 * scale))
            .margin(15 * scale)
            .x_label_area_size(40 * scale)
            .y_label_area_size(50 * scale)
            .build_cartesian_2d(0.0f64..d_x_max * 1.05, 0.0f64..d_y_max * 1.05)?;

        chart.configure_mesh()
            .x_desc("Energy Density \u{3B5} [MeV/fm\u{00B3}]")
            .y_desc("Pressure P [MeV/fm\u{00B3}]")
            .draw()?;

        chart.draw_series(LineSeries::new(clean_data, plot_red.stroke_width(2 * scale)))?
            .label("EoS Model")
            .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], plot_red.stroke_width(2)));

        chart.configure_series_labels().position(SeriesLabelPosition::UpperLeft).draw()?;
    }

    root.present()?;
    Ok(())
}

// ============================================================================
// ARTIST: GERENCIADOR DE MÚLTIPLAS CURVAS (Ideal para comparações)
// ============================================================================

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
    use_logscale: bool, // Controle para escalas logarítmicas
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
            use_logscale: false,
            curves: Vec::new(),
        }
    }

    pub fn with_x_label(mut self, label: &str) -> Self { self.x_label = label.to_string(); self }
    pub fn with_y_label(mut self, label: &str) -> Self { self.y_label = label.to_string(); self }
    pub fn autoscale(mut self) -> Self { self.x_min = None; self.x_max = None; self.y_min = None; self.y_max = None; self }
    pub fn with_x_range(mut self, min: f64, max: f64) -> Self { self.x_min = Some(min); self.x_max = Some(max); self }
    
    /// Habilita escala Log-Log (Ideal para a Equação de Estado)
    pub fn with_log_scale(mut self) -> Self {
        self.use_logscale = true;
        self
    }

    pub fn add_curve(mut self, x: &[f64], y: &[f64], label: &str) -> Self {
        self.curves.push(CurveData { x: x.to_vec(), y: y.to_vec(), label: label.to_string() });
        self
    }

    pub fn plot(&self) -> Result<(), Box<dyn Error>> {
        if self.curves.is_empty() { return Err("Sem curvas.".into()); }

        let path = Path::new(&self.output_path);
        if let Some(parent) = path.parent() { std::fs::create_dir_all(parent)?; }

        let mut d_x_min = f64::MAX; let mut d_x_max = f64::MIN;
        let mut d_y_min = f64::MAX; let mut d_y_max = f64::MIN;

        for curve in &self.curves {
            for (&xi, &yi) in curve.x.iter().zip(curve.y.iter()) {
                if xi.is_finite() && yi.is_finite() {
                    if self.use_logscale && (xi <= 0.0 || yi <= 0.0) { continue; }
                    if !self.use_logscale && self.x_label.contains("Radius") && xi > 50.0 { continue; }
                    
                    if xi < d_x_min { d_x_min = xi; }
                    if xi > d_x_max { d_x_max = xi; }
                    if yi < d_y_min { d_y_min = yi; }
                    if yi > d_y_max { d_y_max = yi; }
                }
            }
        }

        if d_x_min == f64::MAX { d_x_min = 0.1; d_x_max = 1.0; d_y_min = 0.1; d_y_max = 1.0; }

        let x_start = self.x_min.unwrap_or(if self.use_logscale { d_x_min * 0.8 } else if self.x_label.contains("Radius") { d_x_min.max(8.0) } else { 0.0 });
        let x_end = self.x_max.unwrap_or(if self.use_logscale { d_x_max * 1.5 } else { d_x_max * 1.05 });
        let y_start = self.y_min.unwrap_or(if self.use_logscale { d_y_min * 0.8 } else { 0.0 }); 
        let y_end = self.y_max.unwrap_or(if self.use_logscale { d_y_max * 1.5 } else { d_y_max * 1.05 });

        let root = SVGBackend::new(&self.output_path, (800, 600)).into_drawing_area();
        root.fill(&WHITE)?;

        let colors = [&BLUE, &RED, &GREEN, &MAGENTA, &CYAN, &BLACK];

        // A biblioteca plotters requer instâncias separadas de ChartBuilder para Log e Linear
        if self.use_logscale {
            let mut chart = ChartBuilder::on(&root)
                .caption(&self.title, ("sans-serif", 20))
                .margin(20)
                .x_label_area_size(40)
                .y_label_area_size(50)
                .build_cartesian_2d((x_start..x_end).log_scale(), (y_start..y_end).log_scale())?;

            chart.configure_mesh().x_desc(&self.x_label).y_desc(&self.y_label).light_line_style(&WHITE.mix(0.1)).draw()?;

            for (i, curve) in self.curves.iter().enumerate() {
                let color = colors[i % colors.len()];
                let step = (curve.x.len() / 2000).max(1);
                let data = curve.x.iter().zip(curve.y.iter()).enumerate()
                    .filter(|(idx, (x, y))| idx % step == 0 && x.is_finite() && y.is_finite() && *x > &0.0 && *y > &0.0)
                    .map(|(_, (&x, &y))| (x, y));

                chart.draw_series(LineSeries::new(data, color.stroke_width(2)))?
                    .label(&curve.label)
                    .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2)));
            }
            chart.configure_series_labels().position(SeriesLabelPosition::UpperLeft).draw()?;
        } else {
            let mut chart = ChartBuilder::on(&root)
                .caption(&self.title, ("sans-serif", 20))
                .margin(20)
                .x_label_area_size(40)
                .y_label_area_size(50)
                .build_cartesian_2d(x_start..x_end, y_start..y_end)?;

            chart.configure_mesh().x_desc(&self.x_label).y_desc(&self.y_label).light_line_style(&WHITE.mix(0.1)).draw()?;

            for (i, curve) in self.curves.iter().enumerate() {
                let color = colors[i % colors.len()];
                let step = (curve.x.len() / 2000).max(1);
                let data = curve.x.iter().zip(curve.y.iter()).enumerate()
                    .filter(|(idx, (x, y))| idx % step == 0 && x.is_finite() && y.is_finite())
                    .map(|(_, (&x, &y))| (x, y));

                chart.draw_series(LineSeries::new(data, color.stroke_width(2)))?
                    .label(&curve.label)
                    .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2)));
            }
            chart.configure_series_labels().position(SeriesLabelPosition::UpperRight).draw()?;
        }

        root.present()?;
        Ok(())
    }
}