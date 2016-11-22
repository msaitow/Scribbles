
#Pkg.add("Gadfly")
using Gadfly
#Pkg.add("RDatasets")
using RDatasets
#Pkg.add("Cairo")
using Cairo

semi = readtable("semiquinone.csv")

thistheme = Theme(
    panel_stroke=colorant"black",
    major_label_font_size=18pt,
    minor_label_font_size=14pt,
    #key_title_font_size=14pt,
    key_label_font_size=15pt,
    key_position = :right,
    colorkey_swatch_shape = :circle
)

Gadfly.push_theme(thistheme)
Gadfly.set_default_plot_size(24cm, 14cm)

myLineWidth  = 2.0px
myPointWidth = 6.0px

res = plot(semi,
           
           layer(x=:x, y=:DLPNO_CCSD_strong_pairs_, Geom.point, Geom.line,
                 Theme(default_color=colorant"orange",      line_width=myLineWidth, default_point_size=myPointWidth) ),
           
           layer(x=:x, y=:DLPNO_CCSD_total_pairs_,  Geom.point, Geom.line,
                 Theme(default_color=colorant"red",         line_width=myLineWidth, default_point_size=myPointWidth) ),
           
           layer(x=:x, y=:LPNO_CCSD_strong_pairs_,  Geom.point, Geom.line,
                 Theme(default_color=colorant"deepskyblue", line_width=myLineWidth, default_point_size=myPointWidth) ),
           
           layer(x=:x, y=:LPNO_CCSD_total_pairs_,   Geom.point, Geom.line,
                 Theme(default_color=colorant"blue",        line_width=myLineWidth, default_point_size=myPointWidth) ),
           
           Scale.x_log10(), Guide.ylabel("CCSD Correlation Energy Recovered %"), Guide.xlabel("TCutPNO"),
           Guide.manual_color_key("", ["DLPNO-CCSD (strong)", "DLPNO-CCSD (total)", "LPNO-CCSD (strong)", "LPNO-CCSD (total)"], ["orange", "red", "deepskyblue", "blue"]),           
           Coord.Cartesian(ymin=99, ymax=101, xflip=true),
           Stat.xticks(ticks=[-3, -4, -5, -6, -7, -8, -9, -10]),
           Stat.yticks(ticks=[99.0, 99.5, 100.0, 100.5, 101.0]),           
           Guide.xticks(orientation=:horizontal),
           Guide.yticks(ticks=[99.9])           
           #Stat.ygrid(ticks=[99.0, 99.5, 99.6, 99.7, 99.8, 99.9, 100.0, 100.5, 101.0]), 
  )

draw(PDF("semiq.pdf", 25cm, 15cm), res)
