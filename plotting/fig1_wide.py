#65;5802;1c
import svgutils.transform as sg
import sys


def generate_figure(input_julia_path, input_data_path, output_path):
    fig0 = sg.SVGFigure( "210mm", "297mm")
    fig0.root.set("viewBox", "0 0 210 297")
    
    MPL_SCALE = 0.1
    JL_SCALE = 0.08

    def get_file(fn, scale=MPL_SCALE, pos=None):
        fig = sg.fromfile(fn)
        plot = fig.getroot()
        plot.scale_xy(scale, scale)
        if pos is not None:
            plot.moveto(*pos)
        return plot

    def get_mpl_base(fn, scale=MPL_SCALE, pos=None):
        fig = sg.fromfile(fn)
        plot = fig.getroot()
        child = plot[0][0]


        #plot[0] = plot[0][0]
        
        plot.scale_xy(scale, scale)
        if pos is not None:
            plot.moveto(*pos)
        return plot
        
        
    chris_fig = get_file('../extra_panels/HEI10 figure 1.14.svg', scale=1.0)

    fig_panels = [ get_mpl_base(input_data_path+'/rel_length_tot_compare.svg', pos=(14,123)),
                   get_mpl_base(input_data_path+'/new_rel_length.svg', pos=(65,123)),
                   get_mpl_base(input_data_path+'/new_rel_length2.svg', pos=(111,123)),
                   get_mpl_base(input_data_path+'/single_double_triple_peak_new_intensities_stacked.svg', pos=(155,123))]
    
    gpage = sg.GroupElement([chris_fig]+fig_panels)

    fig0.append([gpage])

    # save generated SVG files
    fig0.save(output_path+"/fig1_final.svg")

generate_figure(sys.argv[1], sys.argv[2], sys.argv[3])
