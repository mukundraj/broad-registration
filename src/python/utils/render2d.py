"""
2D rendering helper functions to draw using Cairo

Created by Mukund on 2022-03-22
"""

import cairo
import math

def draw_disk(surface, pts, radius, colors):
    if (len(colors) == 1):
        color = colors[0]
        for i, pt in enumerate(pts):
            x, y, radius = (pt[0], pt[1], radius)
            ctx = cairo.Context(surface)
            ctx.set_line_width(1)
            ctx.arc(x, y, radius, 0, 1.0 * math.pi)
            ctx.set_source_rgb(color[0], color[1], color[2])
            ctx.fill_preserve()
            # ctx.set_source_rgb(1, 1, 1)
            ctx.set_source_rgb(color[0], color[1], color[2])
            ctx.stroke()
    else:
        for i, pt in enumerate(pts):
            x, y, radius = (pt[0], pt[1], radius)
            ctx = cairo.Context(surface)
            ctx.set_line_width(1)
            ctx.arc(x, y, radius, 0, 1.0 * math.pi)
            ctx.set_source_rgb(colors[i][0], colors[i][1], colors[i][2])
            ctx.fill_preserve()
            # ctx.set_source_rgb(1, 1, 1)
            ctx.set_source_rgb(colors[i][0], colors[i][1], colors[i][2])
            ctx.stroke()


def draw_image(pts, width, height, op_file, colors =[[0.5, 0.5, 0.5]], radius=2):

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
    # ctx = cairo.Context(surface)
    # colors =[[0.5, 0.5, 0.5]]
    # radius = 2
    draw_disk(surface, pts, radius, colors)
    surface.write_to_png(op_file)  # Output to PNG
