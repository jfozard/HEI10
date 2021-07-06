

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



    def make_row(data):
        panels = []
        for i, (f, l, c) in enumerate(data):
            panels.append(get_file(f, pos=(i*XS, 0)))
            panels.append(sg.TextElement(i*XS-10, 0, l, size=20, weight="bold"))
            panels.append(sg.TextElement(i*XS+0.5*XS, 0, c, size=24, anchor="middle"))
            
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


    row = [ ( input_data_path+'/violin_number_length.svg', 'a', 'Experimental data'),
            ( input_julia_path+'/new_end_nco_vs_length.svg', 'b', 'Simulation output') ]
    
    g  = make_row(row)
    g.moveto(0,70)



    gpage = sg.GroupElement([g])#, gh])
    gpage.moveto(30,30)


    fig0.append([gpage])



    # save generated SVG files
    fig0.save(output_path+"/CO_length.svg")

generate_figure(*sys.argv[1:])
