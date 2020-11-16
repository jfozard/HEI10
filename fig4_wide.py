

# Underexpressor line

import svgutils.transform as sg
import base64
import sys


def generate_figure(input_julia_path, input_data_path, output_path):
    fig0 = sg.SVGFigure( "210mm", "297mm")


    MPL_SCALE = 0.4
    JL_SCALE = 0.08

    def get_file(fn, scale=MPL_SCALE, pos=None):
        fig = sg.fromfile(fn)
        plot = fig.getroot()
        plot.scale_xy(scale, scale)
        if pos is not None:
            plot.moveto(*pos)
        return plot


    YS = 170
    XS = 185



    def make_row(data, label=True):

        panels = []
        if label:
            for i,f in enumerate(data[0]):
                panels.append(get_file(f, pos=(i*XS, 0)))
            if data[1]:
                panels.append(sg.TextElement(-10, 10, data[1], size=20, weight="bold"))
            if len(data)>2 and data[2]:
                panels.append(sg.TextElement(XS, 0, data[2], size=12, anchor="middle"))
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


    image_panel = sg.ImageElement(open("extra_panels/HEI10 fig 4A.png", "rb"), 1584, 444)
    image_panel.moveto(10,10)
    image_panel.scale_xy(0.48, 0.48)

#    image_label = sg.TextElement(20, 30, "A", size=20, color="white")

    g_image = sg.GroupElement([image_panel])#, image_label])

    data_panels = []


    col1 = [ ((input_data_path+'/ux_all_late_hist_n.svg', input_julia_path+'/ux_new_end_julia_all_hist_n.svg'), 'b', "Bivalent CO number"),
             ((input_data_path+'/ux_all_late_all.svg', input_julia_path+'/ux_new_end_julia_all_hist.svg'), 'c', "All COs"),
             ((input_data_path+'/ux_all_late_1.svg', input_julia_path+'/ux_new_end_julia_hist_1.svg'), 'd', 'Single COs')]


    g  = make_group(col1)
    g.moveto(0,30)


    head1 = sg.TextElement(20, 0, 'Late stage data', size=18)

    head2 = sg.TextElement(XS+10, 0, 'Simulation output', size=18)

    gh = sg.GroupElement([head1, head2])
    gh.moveto(0, 10)


    col2 = [ ((input_data_path+'/ux_all_late_2.svg', input_julia_path+'/ux_new_end_julia_rank_2.svg'), 'e', 'Double COs'),
             ((input_data_path+'/ux_new_rel_peak_int_cell_stacked.svg', input_julia_path+'/ux_new_end_julia_peak_int_stacked.svg'), 'f', 'CO HEI10 focus intensity'),
             ((input_data_path+'/ux_all_spacing.svg', input_julia_path+'/ux_new_end_diff.svg'), 'g', 'CO spacing')]

    g2 = make_group(col2)
    g2.moveto(2*XS, 30)

    head1 = sg.TextElement(20, 0, 'Late stage data', size=18)
    head2 = sg.TextElement(XS+10, 0, 'Simulation output', size=18)

    gh2 = sg.GroupElement([head1, head2])
    gh2.moveto(2*XS, 10)



    gall = sg.GroupElement([g, g2, gh, gh2])
    gall.moveto(30,270)

    gpage = sg.GroupElement([g_image, gall])


    fig0.append([gpage])



    # save generated SVG files
    fig0.save(output_path+"/ux_new_end_fig_final.svg")


generate_figure(sys.argv[1], sys.argv[2], sys.argv[3])
