

# Underexpressor line

import svgutils.transform as sg
import base64
import sys 

def generate_figure(input_julia_path, input_data_path, output_path):
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


    YS = 160*2
    XS = 185*2



    def make_row(data, label=True):

        panels = []
        if label:
            for i,f in enumerate(data[0]):
                panels.append(get_file(f, pos=(i*XS, 0)))
            if data[1]:
                panels.append(sg.TextElement(-10, 10, data[1], size=30, weight="bold"))
            if len(data)>2 and data[2]:
                panels.append(sg.TextElement(XS, 0, data[2], size=20, anchor="middle"))
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
                panels.append(sg.TextElement(label_offset_x, 10, l, size=20, weight="bold"))
        else:
            f=data
            panels.append(get_file(f, pos=(0, 0)))        
        return sg.GroupElement(panels)



    # load matpotlib-generated figures


    col1 = [ ((input_julia_path+'/new_end_julia_all_hist_n_mean.svg', input_julia_path+'/female_julia_all_hist_n_mean.svg'), 'a', "Bivalent CO number"),
             ((input_julia_path+'/new_end_julia_all_hist.svg', input_julia_path+'/female_julia_all_hist.svg'), 'b', "CO Position") ]


    g  = make_group(col1)
    g.moveto(0,70)


    head0 = sg.TextElement(XS, 0, 'Simulation Output', anchor="middle", size=28)

    head1 = sg.TextElement(0.5*XS, 35, 'Male', anchor="middle", size=26)

    head2 = sg.TextElement(1.5*XS, 35, 'Female', anchor="middle", size=26)

    gh = sg.GroupElement([head0, head1, head2])
    gh.moveto(0, 10)


    gpage = sg.GroupElement([g, gh])
    gpage.moveto(30,30)


    fig0.append([gpage])



    # save generated SVG files
    fig0.save(output_path+"/female_male.svg")

generate_figure(*sys.argv[1:])
