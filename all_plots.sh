
base="output/"
julia_dir=${base}"julia_plots/"
kymo_dir=${base}"kymo/"
data_dir="data_output/"
output_dir=${base}"figures/"

mkdir -p ${julia_dir}
mkdir -p ${output_dir}


python julia_plots.py ${base} ${julia_dir}  
python data_analysis_plots.py
python julia_ordered_plots.py ${base} ${julia_dir}
python data_ordered_plots.py

for s in "" no_end_ exp_; do
    u=${s}kymo
    for v in {1..3}; do
	python kymo.py ${kymo_dir}/${u}${v}.dat ${julia_dir}/${u}${v}.svg
    done
    python kymo_tc.py ${kymo_dir}/${u}1.dat ${julia_dir}/${s}timecourse.svg
done


python fig2_wide.py ${julia_dir} ${data_dir} ${output_dir}
python fig3_wide.py ${julia_dir} ${data_dir} ${output_dir}
python fig4_wide.py ${julia_dir} ${data_dir} ${output_dir}
python female_male.py ${julia_dir} ${data_dir} ${output_dir}
python fig_S25.py ${julia_dir} ${data_dir} ${output_dir}
python centromere_fig.py ${julia_dir} ${data_dir} ${output_dir}

inkscape ${output_dir}/new_end_fig_final.svg -o ${output_dir}/fig2.png -d 300 -D
inkscape ${output_dir}/ox_new_end_fig_final.svg -o ${output_dir}/fig3.png -d 300 -D
inkscape ${output_dir}/ux_new_end_fig_final.svg -o ${output_dir}/fig4.png -d 300 -D
inkscape ${output_dir}/no_end_fig_final.svg -o ${output_dir}/figS2.png -d 300 -D
inkscape ${output_dir}/female_male.svg -o ${output_dir}/figS3.png -d 300 -D
inkscape ${output_dir}/centromere_fig.svg -o ${output_dir}/figS4.png -d 300 -D
inkscape extra_panels/figS5.svg -o ${output_dir}/figS5.png -d 300 -D