echo "Empiezan las simulaciones"

matlab -nodisplay -r "publish('k5_D15.m')" &
matlab -nodisplay -r "publish('k10_D15.m')" &
matlab -nodisplay -r "publish('k15_D15.m')" &
matlab -nodisplay -r "publish('k20_D15.m')" &
matlab -nodisplay -r "publish('k25_D15.m')" &
matlab -nodisplay -r "publish('k30_D15.m')" &