

# Code to assemble supplementary figure 4, incorporating centromeric regions

import sys
import svgutils.transform as sg
import base64
import sys

def plot_figure(input_julia_path, input_data_path, output_path):
    fig0 = sg.SVGFigure( "210mm", "297mm")


    MPL_SCALE = 0.4*2
    JL_SCALE = 0.08

    def get_file(fn, scale=MPL_SCALE, pos=None):
        fig = sg.fromfile(fn)
        plot = fig.getroot()
        plot.scale_xy(scale, scale)
        if pos is not None:
            plot.moveto(*pos)
        return plot


    YS = 170*2
    XS = 185*2



    def make_row(data, label=True):

        panels = []
        if label:
            for i,f in enumerate(data[0]):
                panels.append(get_file(f, pos=(i*XS, 0)))
            if data[1]:
                panels.append(sg.TextElement(-10, 10, data[1], size=30, weight="bold"))
            if len(data)>2 and data[2]:
                panels.append(sg.TextElement(XS, 0, data[2], size=24, anchor="middle"))
        else:
            for i,f in enumerate(data):
                panels.append(get_file(f, pos=(i*XS, 0)))        
        return sg.GroupElement(panels)


    def make_group(data, label=True):
        panels = []
        for i,r in enumerate(data):
            g = make_row(r, label)
            g.moveto(0, i*YS)
            panels.append(g)

        return sg.GroupElement(panels)

    def make_single(data, label=True, label_offset_x=-10):
        panels = []
        if label:
            f, l  = data
            panels.append(get_file(f, pos=(0, 0)))
            if l:
                panels.append(sg.TextElement(label_offset_x, 10, l, size=30, weight="bold"))
        else:
            f=data
            panels.append(get_file(f, pos=(0, 0)))        
        return sg.GroupElement(panels)



    # load matpotlib-generated figures


    col1 = [ ((input_data_path+'/ordered_short_n.svg', input_julia_path+'/julia_ordered_short_n.svg'), 'a', "Bivalent CO number"),
             ((input_data_path+'/ordered_short.svg', input_julia_path+'/julia_ordered_short.svg'), 'b', "CO Position") ]


    g  = make_group(col1)
    g.moveto(0,50)



    u1 = make_single((input_data_path+'/short_peak_pos.svg','c'))
    u2 = make_single((input_julia_path+'/ri_plot.svg','d'))
    u2.moveto(XS, 0)
    g2 = sg.GroupElement([u1, u2])
    g2.moveto(0, 2*YS+50)

    head1 = sg.TextElement(0.5*XS, 0, 'Short SC Data', anchor="middle", size=24)

    head2 = sg.TextElement(1.5*XS, 0, 'Simulated Output', anchor="middle", size=24)

    head3 = sg.TextElement(0.5*XS, 2*YS+20, 'DAPI peak position', anchor="middle", size=24)
    head4 = sg.TextElement(1.5*XS, 2*YS+20, 'RI distribution', anchor="middle", size=24)


    gh = sg.GroupElement([head1, head2, head3, head4])
    gh.moveto(0, 10)


    gpage = sg.GroupElement([g, g2, gh])
    gpage.moveto(30,30)


    fig0.append([gpage])



    # save generated SVG files
    fig0.save(output_path+"/centromere_fig.svg")


plot_figure(sys.argv[1], sys.argv[2], sys.argv[3])
