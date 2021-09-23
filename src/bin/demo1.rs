use mui::prelude::*;


fn main() {
    use std::sync::Arc;

    let mut app = App::default();

    let shader = Arc::new(MakeShader2d::make(app.gpu()));
    let white_square = Arc::new(Texture::from_rgba32(app.gpu(), &[u32::MAX], 1, 1));

    let mut draw_op = DrawOp::new(shader);
    draw_op.bind_uniforms::<Sampler2d>(None, Some(&[white_square.clone()]));

    app.insert(Pen2d::new(draw_op));

    app.update(|_app| {

    });

    app.render(|app| {
        let draw_op = {
            let vp_size = app.get_viewport_size().into();
            let mut pen = app.get_mut::<Pen2d>().unwrap();
            pen.viewport_size = vp_size;

            pen.clear();
            pen.fill_rect(Rect::from_xywh(200.0, 150.0, 400.0, 300.0), Color::RED);
            
            pen.bind_to_draw_op();
            pen.draw_op.clone()
        };

        app.backend.renderer.push_draw_op(draw_op);
    });

    run(app);
}
