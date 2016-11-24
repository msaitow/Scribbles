
#Pkg.add("Gadfly")
using Gadfly
#Pkg.add("RDatasets")
using RDatasets
#Pkg.add("Cairo")
using Cairo

semi = readtable("CH2.latest.csv")

thistheme = Theme(
    minor_label_font="AvenirNext-Bold",
    major_label_font="AvenirNext-Bold",
    key_label_font="AvenirNext-Bold",
    
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

xticks = collect(1:1:10)
# for editing
yticks = collect(100.00:0.01:100.05)
# for usual purpose
#yticks = collect(99.0:0.5:101)

res = plot(semi,

           layer(x=:number, y=:scalepair_1_0,  Geom.point, Geom.line,
                 Theme(default_color=colorant"blue",       line_width=myLineWidth, default_point_size=myPointWidth) ),

           layer(x=:number, y=:scalepair_0_5,  Geom.point, Geom.line,
                 Theme(default_color=colorant"orange",     line_width=myLineWidth, default_point_size=myPointWidth) ),
           
           layer(x=:number, y=:scalepair_0_33, Geom.point, Geom.line,
                 Theme(default_color=colorant"red",         line_width=myLineWidth, default_point_size=myPointWidth) ),
                                 
           layer(x=:number, y=:scalepair_0_1,   Geom.point, Geom.line,
                 Theme(default_color=colorant"deepskyblue", line_width=myLineWidth, default_point_size=myPointWidth) ),
           
           Guide.ylabel("CCSD Correlation Energy Recovered %"), Guide.xlabel("Number of CH2 Monomers"),
           Guide.manual_color_key("", ["TScaleTCutPairs=1.00", "TScaleTCutPairs=0.5", "TScaleTCutPairs=0.33", "TScaleTCutPairs=0.1"], ["blue", "orange", "red", "deepskyblue"]),           
           #Coord.Cartesian(ymin=99, ymax=101, xflip=true),
           Stat.xticks(ticks=xticks),
           Stat.yticks(ticks=yticks),           
           Guide.xticks(orientation=:horizontal)
           #Guide.yticks(ticks=[99.9])           
           #Stat.ygrid(ticks=[99.0, 99.5, 99.6, 99.7, 99.8, 99.9, 100.0, 100.5, 101.0]), 
  )

draw(PDF("ch2.pdf", 25cm, 15cm), res)
