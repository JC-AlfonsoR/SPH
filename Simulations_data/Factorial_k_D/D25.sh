echo "Empiezan las simulaciones"

matlab -nodisplay -r "publish('k5_D25.m')" &
matlab -nodisplay -r "publish('k10_D25.m')" &
matlab -nodisplay -r "publish('k15_D25.m')" &
matlab -nodisplay -r "publish('k20_D25.m')" &
matlab -nodisplay -r "publish('k25_D25.m')" &
matlab -nodisplay -r "publish('k30_D25.m')" &