echo "Empieza Simulacion" >> r10.sh
matlab -nodisplay -r "publish('Factorial_k2.m')" &
echo "k2 --- ok" >> r10.sh
matlab -nodisplay -r "publish('Factorial_k4.m')" &
echo "k4 --- ok" >> r10.sh
matlab -nodisplay -r "publish('Factorial_k6.m')" &
echo "k6 --- ok" >> r10.sh
matlab -nodisplay -r "publish('Factorial_k8.m')" &
echo "k8 --- ok" >> r10.sh
matlab -nodisplay -r "publish('Factorial_k10.m')" &
echo "k10 --- ok" >> r10.sh