#!/bin/bash

#SBATCH --job-name=Davi_$SLURM_JOB_ID       # Nome do job
#SBATCH --time=180-00:00:00                 # Tempo máximo de execução (D-hh:mm:ss)
#SBATCH --nodes=1                           # Número de nós
#SBATCH --ntasks-per-node=1                 # Número de tarefas por nó
#SBATCH --cpus-per-task=1                   # Número de CPUs por tarefa
#SBATCH --output=meujob-%j.out              # Arquivo de saída do slurm
#SBATCH --mail-type=ALL                     # Tipo de notificação por email
#SBATCH --mail-user=davi.icex@gmail.com     # Endereço de email para notificação


# Parâmetros da simulação
theta=0.
nx=12
contorno=1
tipo=2
copias=3
n_mc=100000
n_temp=50
ti=20.
tf=1.
single=10
kago=5
triang=0

# Cria o arquivo com os parâmetros da simulação
parametros="parametros_${theta}.dat"
echo "! Parâmetros da rede a ser simulada !" > "$parametros"
echo "${theta},        !dTheta: Rotação dos spins" >> "$parametros"
echo "${nx},           !Nx: Repetições da célular unitária em x" >> "$parametros"
echo "${nx},           !Ny: Repetições da célular unitária em y" >> "$parametros"
echo "${contorno},     !Contorno: 1 -> com contorno, 2 -> sem contorno" >> "$parametros"
echo "${tipo},         !Troca: 0 -> sem troca, 1 -> com troca, 2 -> D = 1.0" >> "$parametros"
echo "${copias},       !nc: Número de cópias" >> "$parametros"
echo "${n_mc},         !N_mc: Número de passos de Monte Carlo" >> "$parametros"
echo "${n_temp},       !N_temp: Números de pontos de Temperaturas" >> "$parametros"
echo "${ti},           !Ti: Temperatura inicial" >> "$parametros"
echo "${tf},           !Tf: Temperatura final" >> "$parametros"
echo "${single},       !N_single: Passos de single-spin-flip" >> "$parametros"
echo "${kago},         !N_kago: Passos de loops kagome" >> "$parametros"
echo "${triang},       !N_tri: Passos de loops triangular" >> "$parametros"


# Cria um novo diretório para este job com base na id do job
DIR="theta_${theta}_job_${SLURM_JOB_ID}"
mkdir $DIR
cd $DIR

# Copia o código Fortran para o novo diretório
cp ../MC_in.f90 .
cp ../"$parametros" ./"parametros_rede.dat"

# Compila o código usando gfortran
gfortran -o MC MC_in.f90

# Executa o programa
./MC

# Volta para o diretório inicial
cd ..

# Compacta a pasta gerada
tar -cvzf output_theta_${theta}_${SLURM_JOB_ID}.tar.gz $DIR

# Exclui o diretório do job
rm -rf $DIR
