

## MUESTREADOR PERSONALIZADO METROPOLIS-HASTINGS
## muestreador para actualizar conjuntamente y.un[1:M,j] de manera que ponemos
## en cada paso del muestreo la condicion de que sumen n[j]
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Definimos los componentes a usar
    nnidd<-control$nnidd
    j<-control$j
    M<-control$M
    calcNodes <- model$getDependencies(target)
  },

run = function() {
  lam.curr <- model$lam[1:M,j] # Conteos esperados de individuos por trampa

  #Muestrea y[1:M,j] reasignando n[j] usando el condicional completo
  switch.probs <- lam.curr[1:M]/sum(lam.curr[1:M])

  #propone nuevas identificaciones para nnid[j,k]
  y.latent.curr <- model$y.full[1:M,j]- model$y.obs[1:M,j]
  y.latent.prop <- rmulti(1, nnidd, switch.probs[1:M])
  model$y.full[1:M,j] <<-  model$y.obs[1:M,j] + y.latent.prop

  # modelo inicial logProb
  model_lp_initial <- model$getLogProb(calcNodes)

  # modelo propuesto logProb
  model_lp_proposed <- model$calculate(calcNodes)

  # Relación log-Metropolis-Hastings
  log_MH_ratio <- (model_lp_proposed + dmulti(y.latent.curr, nnidd, switch.probs, log=TRUE)) -
                  (model_lp_initial + dmulti(y.latent.prop,  nnidd, switch.probs, log=TRUE))

  # Paso Metrópolis-Hastings
  accept <- decide(log_MH_ratio)
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list( reset = function () {} )
)
