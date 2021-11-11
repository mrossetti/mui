use mui::prelude::*;


fn main() {
    use std::sync::Arc;

    let mut app = App::default();
    
    let font_atlas = {
        if std::path::Path::new("local/font_atlas.ron").exists() {
            println!("Attempting to load FontAtlas...");
            FontAtlas::load(app.gpu(), "local/font_atlas.png")
        } else {
            let proggy_clean = ab_glyph_support::load_font("res/fonts/ProggyClean.ttf");
            let font_atlas = ab_glyph_support::make_font_atlas(app.gpu(), &[(proggy_clean, Default::default(), CharacterSet::ASCII)]);
            font_atlas.save("local/font_atlas.png");
            font_atlas
        }
    };

    let shader = Arc::new(MakeShader2dR::make(app.gpu()));
    let mut draw_op = DrawOp::new(shader);
    draw_op.bind_uniforms::<Sampler2d>(&[], &[Arc::clone(&font_atlas.texture)]);
    let stream = Stream::<Vertex4fx2>::new(draw_op);

    pub struct Pen {
        pub stream: Stream<Vertex4fx2>,
        pub font_atlas: FontAtlas,
    }

    app.insert(Pen { stream, font_atlas });

    // main loop callbacks
    app.update(|_app| {

    });

    app.render(|app| {
        let draw_op = {
            let vp_size = app.get_viewport_size().into();
            let pen = app.get_mut::<Pen>().unwrap();
            pen.stream.viewport_size = vp_size;

            pen.stream.clear();
            v2d::draw_text(&mut pen.stream, &mut pen.font_atlas, "Hello, world!", 
            DrawTextArgs {
                px_dest: vec2(20.0, 20.0),
                scale: vec2(2.0, 2.0),
                underline: Underline::THICK4,
                ..Default::default()
            });
            
            pen.stream.update_draw_op_buffers();
            pen.stream.draw_op.clone()
        };

        app.backend.renderer.push_draw_op(draw_op);
    });

    run(app);
}
