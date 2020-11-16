#65;5802;1c
import svgutils.transform as sg
import sys

def plot_prefix(prefix, julia_input_path, data_input_path, output_path):
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


    YS = 150
    YS = 150
    XS = 200

    #         (('rel_peak_int_cell_stacked.svg', 'ux_new_end_julia_peak_int_stacked.svg'), 'H'),

    col1 = [ ((data_input_path+'/all_late_hist_n.svg', julia_input_path+'/'+prefix+'julia_all_hist_n.svg'), 'a', "Bivalent CO number"),
             ((data_input_path+'/all_late_all.svg', julia_input_path+'/'+prefix+'julia_all_hist.svg'), 'b', "All COs"),
             ((data_input_path+'/all_late_1.svg', julia_input_path+'/'+prefix+'julia_hist_1.svg'), 'c', "Single COs"),
             ((data_input_path+'/all_late_2.svg', julia_input_path+'/'+prefix+'julia_rank_2.svg'), 'd', "Double COs"),
             ((data_input_path+'/all_late_3.svg', julia_input_path+'/'+prefix+'julia_rank_3.svg'), 'e', "Triple COs"),
             ((data_input_path+'/all_spacing.svg', julia_input_path+'/'+prefix+'diff.svg'), 'f', "CO spacing") ]



    #         

    def make_row(data, label=True, ls=12):

        panels = []
        if label:
            for i,f in enumerate(data[0]):
                panels.append(get_file(f, pos=(i*XS, 0)))
            if data[1]:
                panels.append(sg.TextElement(-10, 10, data[1], size=20, weight="bold"))
            if len(data)>2 and data[2]:
                panels.append(sg.TextElement(XS, 7, data[2], size=12, anchor="middle"))
        else:
            for i,f in enumerate(data):
                panels.append(get_file(f, pos=(i*XS, 0)))        
        return sg.GroupElement(panels)


    def make_group(data, label=True, ls=12):
        panels = []
        for i,r in enumerate(data):
            g = make_row(r, label, ls=ls)
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

    g = make_group(col1)

    g.moveto(0,30)


    head1 = sg.TextElement(25, 0, 'Late stage data', size=18)
    head2 = sg.TextElement(XS+0+10, 0, 'Simulation output', size=18)
    gh = sg.GroupElement([head1, head2])
    gh.moveto(0, 20)


    OT = 45

    # add text labels



    block = [ ((julia_input_path+'/'+prefix+'kymo1.svg', julia_input_path+'/'+prefix+'kymo2.svg'), 'h') , (( julia_input_path+'/'+prefix+'kymo3.svg', julia_input_path+'/'+prefix+'timecourse.svg'),'') ]


    g2 = make_group(block)

    #

    g2.moveto(2*XS, 1.7*YS+30)


    txt2 = [ sg.TextElement(30, OT+i*YS, t, size=16, anchor="middle") for i,t in enumerate(["Simulation time course"])]


    gt2 = sg.GroupElement(txt2)
    gt2.moveto(2*XS+XS-35, 1.7*YS-10)



    plot_hs = make_single((julia_input_path+'/'+prefix+"julia_hs.svg","g"), label_offset_x=-40)
    head_hs = sg.TextElement(20, 0, 'CO Homeostasis', size=16)

    group_hs = sg.GroupElement([plot_hs, head_hs])

    group_hs.moveto(2.5*XS, 0.7*YS)


    col3 = [ ((data_input_path+'/new_centre.svg', julia_input_path+'/'+prefix+'thinned_centre.svg'), 'i', "Double CO HEI10 focus intensity vs position"),
             ((data_input_path+'/new_rel_peak_int_cell_stacked.svg', julia_input_path+'/'+prefix+'julia_peak_int_stacked.svg'), 'j', "CO HEI10 focus intensity") ]

    head1 = sg.TextElement(25, 0, 'Late stage data', size=18)
    head2 = sg.TextElement(XS+0+10, 0, 'Simulation output', size=18)
    gh2 = sg.GroupElement([head1, head2])
    gh2.moveto(2*XS, 4*YS+0)

    g3 = make_group(col3)
    g3.moveto(2*XS, 4*YS+30)

    g3[0][3].moveto(0, -10)
    print(g3[0][3])

    g3[0][2].moveto(0, -20)


    gall = sg.GroupElement([g, g2, g3, gh, gt2, gh2, g3, group_hs])

    gall.scale_xy(0.9, 0.9)
    gall.moveto(30,0)

    gpage = sg.GroupElement([gall])


    fig0.append([gpage])



    # save generated SVG files
    fig0.save(output_path+'/'+prefix+"fig_final.svg")

def generate_figures(julia_input_path, data_input_path, output_path):
    plot_prefix('no_end_', julia_input_path, data_input_path, output_path)
    plot_prefix('exp_', julia_input_path, data_input_path, output_path)

generate_figures(sys.argv[1], sys.argv[2], sys.argv[3])
