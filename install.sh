# https://github.com/YosefLab/Cassiopeia
export GUROBI_HOME="/users/PAS1571/wangcankun100/Cassiopeia/gurobi911/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

# Start interactive job at osc
salloc --nodes=1 --ntasks=8 --mem=64GB --account PCON0022 --time=12:00:00 srun --pty /bin/bash