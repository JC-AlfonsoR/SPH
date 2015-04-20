echo "Empieza Simulacion" >> r20.sh
matlab -nodisplay -r "publish('Factorial_k12.m')" &
echo "k12 --- ok" >> r20.sh
matlab -nodisplay -r "publish('Factorial_k14.m')" &
echo "k14 --- ok" >> r20.sh
matlab -nodisplay -r "publish('Factorial_k16.m')" &
echo "k16 --- ok" >> r20.sh
matlab -nodisplay -r "publish('Factorial_k18.m')" &
echo "k18 --- ok" >> r20.sh