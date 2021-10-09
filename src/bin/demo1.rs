use mui::prelude::*;


fn main() {
    use std::sync::Arc;

    let mut app = App::default();

    let shader = Arc::new(MakeShaderText2d::make(app.gpu()));
    // let white_square = Arc::new(Texture::from_rgba32(app.gpu(), &[u32::MAX], 1, 1));

    #[allow(unused_mut)]
    let mut draw_op = DrawOp::new(shader);
    // draw_op.bind_uniforms::<Sampler2d>(None, Some(&[white_square.clone()]));

    let font = Font::load(app.gpu(), "res/fonts/ProggyClean.ttf", true);

    save_bitmap("local/font_atlas.png", &font.atlas.texture_data, font.atlas.texture_size.x as _, font.atlas.texture_size.y as _);

    let mut stream = Stream2d::new(draw_op);
    {
        let font_atlas_texture = font.atlas.texture.clone();
        stream.bind_font::<Sampler2d>(Box::new(font), None, &[font_atlas_texture]);
    }

    app.insert(stream);

    // main loop callbacks
    app.update(|_app| {

    });

    app.render(|app| {
        let draw_op = {
            let vp_size = app.get_viewport_size().into();
            let mut stream = app.get_mut::<Stream2d>().unwrap();
            stream.viewport_size = vp_size;

            stream.clear();
            stream.push_text("Hello, world!", &Default::default());
            
            stream.bind_to_draw_op();
            stream.draw_op.clone()
        };

        app.backend.renderer.push_draw_op(draw_op);
    });

    run(app);
}
