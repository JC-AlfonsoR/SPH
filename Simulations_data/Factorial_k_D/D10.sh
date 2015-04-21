echo "Empiezan las simulaciones"

matlab -nodisplay -r "publish(k5_D10.m)" &
matlab -nodisplay -r "publish(k10_D10.m)" &
matlab -nodisplay -r "publish(k15_D10.m)" &
matlab -nodisplay -r "publish(k20_D10.m)" &
matlab -nodisplay -r "publish(k25_D10.m)" &
matlab -nodisplay -r "publish(k30_D10.m)" &