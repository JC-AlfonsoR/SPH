echo "Empiezan las simulaciones"

matlab -nodisplay -r "publish(k5_D5.m)" &
matlab -nodisplay -r "publish(k10_D5.m)" &
matlab -nodisplay -r "publish(k15_D5.m)" &
matlab -nodisplay -r "publish(k20_D5.m)" &
matlab -nodisplay -r "publish(k25_D5.m)" &
matlab -nodisplay -r "publish(k30_D5.m)" &