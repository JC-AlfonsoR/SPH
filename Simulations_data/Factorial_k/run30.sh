echo "Empieza Simulacion" >> r30.sh
matlab -nodisplay -r "publish('Factorial_k20.m')" &
echo "k20 --- ok" >> r30.sh
matlab -nodisplay -r "publish('Factorial_k25.m')" &
echo "k25 --- ok" >> r30.sh
matlab -nodisplay -r "publish('Factorial_k30.m')" &
echo "k30 --- ok" >> r30.sh
