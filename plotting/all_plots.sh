
set -e

base="../output/"
julia_dir=${base}"julia_plots/"
kymo_dir=${base}"kymo/"
data_dir=${base}"/data_output/"
output_dir=${base}"figures/"

mkdir -p ${julia_dir}
mkdir -p ${output_dir}



python julia_plots.py ${base} ${julia_dir}  


python data_analysis_plots.py ${base}"/data_output/"

python SC_length_vs_CO_number.py ${data_dir} ${data_dir}


for s in "" no_end_ exp_ escape_; do
    u=${s}kymo
    for v in {1..3}; do
	python kymo.py ${kymo_dir}/${u}${v}.dat ${julia_dir}/${u}${v}.svg
    done
    python kymo_tc.py ${kymo_dir}/${u}1.dat ${julia_dir}/${s}timecourse.svg
done

python fig1_wide.py ${julia_dir} ${data_dir} ${output_dir}
python fig2_wide.py ${julia_dir} ${data_dir} ${output_dir}
python fig3_wide.py ${julia_dir} ${data_dir} ${output_dir}
python fig4_wide.py ${julia_dir} ${data_dir} ${output_dir}
python female_male.py ${julia_dir} ${data_dir} ${output_dir}
python fig_S25.py ${julia_dir} ${data_dir} ${output_dir}

python CO_length.py ${julia_dir} ${data_dir} ${output_dir}


inkscape ${output_dir}/fig1_final.svg -o ${output_dir}/fig1.png -d 300 -D
inkscape ${output_dir}/new_end_fig_final.svg -o ${output_dir}/fig2.png -d 300 -D
inkscape ${output_dir}/ox_new_end_fig_final.svg -o ${output_dir}/fig3.png -d 300 -D
inkscape ${output_dir}/ux_new_end_fig_final.svg -o ${output_dir}/fig4.png -d 300 -D
inkscape ${output_dir}/no_end_fig_final.svg -o ${output_dir}/figS3.png -d 300 -D
inkscape ${output_dir}/female_male.svg -o ${output_dir}/figS4.png -d 300 -D
inkscape ${output_dir}/CO_length.svg -o ${output_dir}/figS5.png -d 300 -D
inkscape ${output_dir}/escape_fig_final.svg -o ${output_dir}/figS6.png -d 300 -D
