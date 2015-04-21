echo "Empiezan las simulaciones"

matlab -nodisplay -r "publish(k5_D30.m)" &
matlab -nodisplay -r "publish(k10_D30.m)" &
matlab -nodisplay -r "publish(k15_D30.m)" &
matlab -nodisplay -r "publish(k20_D30.m)" &
matlab -nodisplay -r "publish(k25_D30.m)" &
matlab -nodisplay -r "publish(k30_D30.m)" &