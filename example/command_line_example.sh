for EQ in 2.50 2.70 2.90 3.10 3.30 3.50 3.70 3.90 4.10 4.30 4.50 4.70 4.90 5.10
do
    python ../format_input_data.py -i data/com_distance_vs_time_files/Nafc500_pmf_pullx_${EQ}_displacment_reduced.dat -f 500 -r $EQ -o data/com_distance_vs_time_files/Nafc500_${EQ}.ui_dat
done
python ../umbrella_integration.py -t 298.15 -i data/com_distance_vs_time_files/Nafc500*.ui_dat -b 0.2 -r right -m 2.5 5.1 -o example.pmf -nb 15 -ph position_histograms.png
