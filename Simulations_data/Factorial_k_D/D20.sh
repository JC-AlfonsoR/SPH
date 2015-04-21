echo "Empiezan las simulaciones"

matlab -nodisplay -r "publish(k5_D20.m)" &
matlab -nodisplay -r "publish(k10_D20.m)" &
matlab -nodisplay -r "publish(k15_D20.m)" &
matlab -nodisplay -r "publish(k20_D20.m)" &
matlab -nodisplay -r "publish(k25_D20.m)" &
matlab -nodisplay -r "publish(k30_D20.m)" &