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


/// Estrutura para armazenar os dados individuais de cada curva
#[derive(Clone)]
pub struct CurveData {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub label: String,
}

/// O gerenciador de gráficos (Plotter)
pub struct Artist {
    output_path: String,
    title: String,
    x_label: String,
    y_label: String,
    x_min: Option<f64>,
    x_max: Option<f64>,
    curves: Vec<CurveData>,
}

impl Artist {
    /// Inicializa um novo Plotter vazio
    pub fn new(output_path: &str, title: &str) -> Self {
        Self {
            output_path: output_path.to_string(),
            title: title.to_string(),
            x_label: "Radius [km]".to_string(),
            y_label: "Mass [M\u{2299}]".to_string(),
            x_min: None,
            x_max: None,
            curves: Vec::new(),
        }
    }

    /// Permite customizar o rótulo do eixo X
    pub fn with_x_label(mut self, label: &str) -> Self {
        self.x_label = label.to_string();
        self
    }

    /// Define manualmente o intervalo do eixo X (ex.: fixar raio 8-14 km)
    pub fn with_x_range(mut self, min: f64, max: f64) -> Self {
        self.x_min = Some(min);
        self.x_max = Some(max);
        self
    }

    /// Permite customizar o rótulo do eixo Y
    pub fn with_y_label(mut self, label: &str) -> Self {
        self.y_label = label.to_string();
        self
    }

    /// Adiciona uma nova curva ao gráfico. Usa o padrão Builder.
    pub fn add_curve(mut self, x: &[f64], y: &[f64], label: &str) -> Self {
        self.curves.push(CurveData {
            x: x.to_vec(),
            y: y.to_vec(),
            label: label.to_string(),
        });
        self
    }

    /// Executa a renderização do gráfico com todas as curvas adicionadas
    pub fn plot(&self) -> Result<(), Box<dyn Error>> {
        if self.curves.is_empty() {
            return Err("Nenhuma curva foi adicionada ao Plotter.".into());
        }

        let path = Path::new(&self.output_path);
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }

        // --- Configuração High DPI ---
        let scale = 1; 
        let width = 800 * scale;
        let height = 600 * scale;

        let root = SVGBackend::new(&self.output_path, (width, height)).into_drawing_area();
        root.fill(&WHITE)?;

        // Encontrar os limites globais dinâmicos para os eixos baseados em todas as curvas
        let mut x_min = self.x_min.unwrap_or(f64::MAX);
        let mut x_max = self.x_max.unwrap_or(f64::MIN);
        let mut y_max = f64::MIN;

        for curve in &self.curves {
            let local_x_min = curve.x.iter().fold(f64::MAX, |a, &b| a.min(b));
            let local_x_max = curve.x.iter().fold(f64::MIN, |a, &b| a.max(b));
            let local_y_max = curve.y.iter().fold(f64::MIN, |a, &b| a.max(b));

            if self.x_min.is_none() && local_x_min < x_min { x_min = local_x_min; }
            if self.x_max.is_none() && local_x_max > x_max { x_max = local_x_max; }
            if local_y_max > y_max { y_max = local_y_max; }
        }

        // Adiciona margens de respiração (5% acima e abaixo) somente quando o eixo é dinâmico
        if self.x_min.is_none() { x_min *= 0.95; }
        if self.x_max.is_none() { x_max *= 1.05; }
        y_max *= 1.05;

        // Se o x_min for muito pequeno, forçamos um limite para focar na estrela (ex: 8.0km)
        if self.x_min.is_none() && x_min < 8.0 && self.x_label.contains("Radius") {
            x_min = 8.0; 
        }

        let mut chart = ChartBuilder::on(&root)
            .caption(&self.title, ("sans-serif", 24 * scale).into_font())
            .margin(15 * scale)
            .x_label_area_size(40 * scale)
            .y_label_area_size(50 * scale)
            .build_cartesian_2d(x_min..x_max, 0.0f64..y_max)?;

        chart
            .configure_mesh()
            .x_desc(&self.x_label)
            .y_desc(&self.y_label)
            .axis_desc_style(("sans-serif", 16 * scale))
            .draw()?;

        // Paleta de cores para diferenciar as curvas
        let colors = [
            &BLUE, &RED, &GREEN, &MAGENTA, &CYAN, &BLACK, &YELLOW
        ];

        // Desenha cada curva registrada
        for (i, curve) in self.curves.iter().enumerate() {
            let color = colors[i % colors.len()];
            let data_iter = curve.x.iter().copied().zip(curve.y.iter().copied());

            chart
                .draw_series(LineSeries::new(
                    data_iter,
                    color.stroke_width(2 * scale),
                ))?
                .label(&curve.label)
                .legend({
                    let c = *color; // Copy color for the closure
                    move |(x, y)| PathElement::new(
                        vec![(x, y), (x + (20 * scale) as i32, y)], 
                        c.stroke_width(2 * scale)
                    )
                });
        }

        chart
            .configure_series_labels()
            .position(SeriesLabelPosition::UpperRight)
            .background_style(&WHITE.mix(0.8))
            .border_style(&BLACK)
            .label_font(("sans-serif", 12 * scale))
            .draw()?;

        root.present()?;
        Ok(())
    }
}