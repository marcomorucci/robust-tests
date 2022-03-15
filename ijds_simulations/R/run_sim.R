library(argparser)
source('simulations.R')

parser <- arg_parser("Simulations for Robust tests for matching")

parser <- add_argument(parser, "simulation", help="Which simulation to run")
parser <- add_argument(parser, "--outpath", help="Output folder", 
                       default="./")
parser <- add_argument(parser, "--dgp", help="Outcome data dgp", 
                       default="gen_linear")
parser <- add_argument(parser, "--samplesize", help="sample sizes",
                       default="50 100 250 500")
parser <- add_argument(parser, "--covariates", help="number of covariates", 
                       default="20")
parser <- add_argument(parser, "--ate", help="treatment effects", 
                       default="0.0 0.2 0.5 0.8 1.2 2")
parser <- add_argument(parser, "--proptreated", help="proportion treated", 
                       default="0.2")
parser <- add_argument(parser, "--var", help="outcome variance", 
                       default="0.1")

argv <- parse_args(parser)

sim_fcn <- get(argv$simulation)
out_path <- argv$outpath
dgp <- get(argv$dgp)
samplesize <- as.numeric(unlist(strsplit(argv$samplesize, " ")))
covariates <- as.numeric(unlist(strsplit(argv$covariates, " ")))
ate <- as.numeric(unlist(strsplit(argv$ate, " ")))
proptreated <- as.numeric(unlist(strsplit(argv$proptreated, " ")))
var <- as.numeric(unlist(strsplit(argv$var, " ")))


# Global options
t_max = 60*5
solver = "cplex"
approximate = 1
solver = list(name = solver, t_max = t_max, approximate = approximate,
              round_cplex = 0, trace_cplex = 0)
cplex_control = list(trace = 0, round = 0, tilim = 60*5, threads=1)
kt = 1
kc = 1
eps = 1e-1
approximate = TRUE
tols = seq(0.01, 1, by=0.01)


out  <- run_sim_slurm(sim_fcn, dgp, samplesize, ate, covariates, proptreated, var)

.rslurm_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

saveRDS(out, file = paste0(out_path, '/results_', .rslurm_id, '.RDS'))




