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
    let r_max = radii_km.iter().fold(0./0., |a: f64, b| a.max(*b)) * 1.05;
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