### Script para SDM multitemporal ###

### abrindo pacotes
library(raster)
library(biomod2)
library(sensitivity)
library(ENMeval)
library(dismo)
library(usdm)
library(rgdal)


### ajustando nosso diretorio de trabalho
# setwd('/home/anderson/Área de Trabalho/sdm_multitemp')


### carregando nossas funcoes
source('./utils/paleoextract.R')
source('./utils/strings2na.R')
source('./utils/dataInstance.R')
source('./utils/uplot.R')
source('./utils/uniche.R')
source('./utils/cleanData.R')
source('./utils/paleobg.R')
source('./utils/ustpv.R')


### definindo diretorios que iremos utilizar
projectFolder = getwd()  # pasta de trabalho do projeto
envFolder = file.path(projectFolder, 'environmental_data', 'climate_data')  # pasta das variaveis ambientais
maxentFolder = file.path(projectFolder,'maxent')  # caminho ate o arquivo do Maxent
AmSulBorders = rgdal::readOGR(
  file.path(projectFolder, 'environmental_data', 'geographical_data', 'borders.shp')
)  # shapefile da Am. do Sul


### selecao (multitemporal) de variaveis
# stpvTable = ustpv(path = envFolder, ageMin = 0, ageMax = 5)  # intervalo entre 1 e 5 mil anos antes do presente

# head(stpvTable) # inspecionando a tabela

# apenas variaveis com todas ou quase todas as linhas com valor = 1
preditorsStpv = c('bio3_stham', 'bio8_stham', 'bio13_stham', 'bio14_stham','bio15_stham','bio18_stham','bio19_stham')

dataSetRaw = read.csv(
    file = file.path(projectFolder, 'occurrence_data', 'dataset.csv'),
  header = TRUE, 
  dec = '.', 
  sep = ','
)  # abrindo

head(dataSetRaw)  # inspecionando visualmente


### ajustes preparatorio nos dados
# transformando strings ao longo do dateset em NA (por exemplo: 
# "gruta da Fazenda Esmeralda" ? transformado em NA)
pts = strings2na(x = dataSetRaw, safeCols = 'Species')
names(pts) <- c("Species", "lon", "lat", "ageMean", "ageMin", "ageMax")


### rodando analise de incerteza
# cols = c("lon", "lat", "ageMean", "ageMin", "ageMax", preditorsStpv) # nomes das colunas do arquivo

# analise.incerteza = uniche(x=pts, cols=cols, dataMaxAge=5, envFolder=envFolder, maxentFolder=maxentFolder)

# uplot(analise.incerteza, AmSulBorders, legend=FALSE) #inspecionando resultados da analise de incerteza
# dev.off()

### pegando dataset
# OBS.: se precisar retirar algum ponto, seguir esse código de exemplo:
#ptsToExclude = c(9, 15, 29, 31, 32)
#occData = analise.incerteza$dataset[ -match(ptsToExclude, analise.incerteza$dataset$ID),]

# occData = analise.incerteza$dataset  # aqui, para o caso de ter rodado o uniche
occData = pts  # aqui, para o caso de NÃO ter rodado o uniche

### gerando instancia de dados para a modelagem
occDataList = dataInstance(x=occData, col_names=c('ageMean', 'ageMin', 'ageMax'))


### rodando SDMs
for(i in 1:length(occDataList)){

  # dados de ocorrencia
  occDataList[[i]]$id = 1:nrow(occDataList[[i]])
  current_occData = paleoextract(x=occDataList[[i]], path=envFolder)
  current_occData = current_occData[complete.cases(current_occData),]
  
  # dados de pseudo-ausencia
  current_bgData = paleobg(x=current_occData, colNames=c('lon','lat','age'), envFolder=envFolder, n=10000) #criando pontos de fundo (i.e., backfround points)
  current_bgData[,c('lon','lat')] = round(current_bgData[,c('lon','lat')], 2) #arredondando de acordo com a sua necessiadade
  current_bgData = current_bgData[!duplicated(current_bgData[,c('lon','lat','age')]), ] #excluindo pontos duplicados (dentro de uma mesma idade)
  current_bgData$id = NA
  current_bgData = current_bgData[names(current_occData)]
  
  # ocorrencias + pseudo-ausencias
  current_occData$type = 1 #presencas = 1
  current_bgData$type = 0 #pseudo-ausencias = 0
  current_dataSet = rbind(current_occData, current_bgData) #dataset consolidado
  current_dataSet = current_dataSet[complete.cases(current_dataSet[, grep(pattern='id', x=names(current_dataSet), invert=TRUE) ]), ] #uma lipezinha basica  : )
  
  # variaveis e parametros locais especificos para o biomod2 (passando para o mundo do biomod2!)
  myRespName = paste('instancia_', i, sep='')
  myResp = current_dataSet$type
  myRespXY = current_dataSet[,c('lon','lat')]
  myExpl = current_dataSet[,grep(pattern="lon|lat|id|age|type", x=names(current_dataSet),invert=TRUE)] 

  # ajuste de dados de entrada para biomod2
  myBiomodData = BIOMOD_FormatingData(resp.var = myResp,
                                      expl.var = myExpl,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

  # parametrizando os modelos
  myBiomodOption = BIOMOD_ModelingOptions(
    MAXENT.Phillips=list(
      path_to_maxent.jar      = maxentFolder,
      linear                  = TRUE,
      quadratic               = TRUE,
      product                 = FALSE,
      threshold               = FALSE,
      hinge                   = FALSE,
      doclamp                 = FALSE,
      maximumiterations       = 1000,
      convergencethreshold    = 1.0E-5,
      threads=4)
  )

  setwd(paste(projectFolder, '/outputs', sep=''))

  # rodando o algoritmo (i.e. SDMs)
  myBiomodModelOut = BIOMOD_Modeling(
    myBiomodData,
    models            = c('MAXENT.Phillips'),
    bm.options        = myBiomodOption,
    nb.rep            = 3,
    data.split.perc   = 75,
    var.import        = 0,
    metric.eval       = c('TSS','ROC'),
    save.output       = FALSE,
    scale.models      = FALSE,
    do.full.models    = FALSE,
    modeling.id       = myRespName
  )
  
    setwd(projectFolder)

  # dados de avaliacao do modelo
  evaluationScores = get_evaluations(myBiomodModelOut)
  
  # salvando dados no HD
  save(myBiomodData, file=paste(projectFolder, '/outputs/', gsub('_', '.', myRespName), '/', gsub('_', '.', myRespName), '.myBiomodData.R', sep=''))
  save(myBiomodOption, file=paste(projectFolder, '/outputs/', gsub('_', '.', myRespName), '/', gsub('_', '.', myRespName),'.myBiomodOption.R',sep=''))
  save(myBiomodModelOut, file=paste(projectFolder, '/outputs/', gsub('_', '.', myRespName), '/', gsub('_', '.', myRespName),'.myBiomodModelOut.R',sep=''))

}


### implementando projecoes
for(i in seq(length(occDataList))){
  for(year in 1:5){  # trabalhando com intervalo entre 1 e 5 mil anos, como exemplo
    
    envVarPaths = list.files( file.path(envFolder, year), full.names=TRUE, pattern='.asc') # vamos projetar para 5  mil anos atras (ultimo maximo glacial)
    predictors = stack(envVarPaths) #abrindo e armazenando as variaveis para dentro do R
    myRespName = paste('instancia_', i, sep='')
    
    # carregando dados
    load(paste(projectFolder, '/outputs/', gsub('_', '.', myRespName), '/', gsub('_', '.', myRespName), '.myBiomodData.R', sep=''))
    load(paste(projectFolder,'/outputs/', gsub('_', '.', myRespName), '/', gsub('_', '.', myRespName),'.myBiomodOption.R',sep=''))
    #load(paste(projectFolder,'/', gsub('_', '.', myRespName), '/', gsub('_', '.', myRespName),'.myBiomodModelOut.R',sep=''))
    
    setwd(paste(projectFolder, '/outputs/projections', sep=''))

    # calibrando um modelo usando todo o conjunto de dados
    myBiomodModelFull = BIOMOD_Modeling(
      myBiomodData,
      models            = c('MAXENT.Phillips'),
      bm.options        = myBiomodOption,
      nb.rep            = 1,
      data.split.perc   = 100,
      var.import        = 0,
      metric.eval       = c('TSS','ROC'),
      save.output       = FALSE,
      scale.models      = FALSE,
      do.full.models    = FALSE,
      modeling.id       = paste(myRespName,'_FULL',sep='')
    )
    
    # rodando algoritmo de projecao
    myBiomodProj = BIOMOD_Projection(
      bm.mod                = myBiomodModelFull,
      new.env               = predictors,
      proj.name             = paste(year,'kyr',sep=''),
      models.chosen         = 'all',
      metric.binary         = c('TSS','ROC'),
      compress              = 'TRUE',
      build.clamping.mask   = 'TRUE',
      output.format         = '.grd'
    )
    
    setwd(projectFolder)

    # raster para as projecoes
    rasterData = get_predictions(myBiomodProj)
    rasterData = mean(rasterData)
    writeRaster(
      x = rasterData,
      filename = paste(projectFolder, '/outputs/projections/', gsub('_', '.', myRespName), '/proj_', year, 'kyr.asc', sep=''),
      overwrite = TRUE
    )
  }
}


### analisando as projecoes
number_of_instances = 1  # ajustar aqui, de acordo com o que tiver sido usado nas etapas anteriores 
number_of_projections = 1  # ajustar aqui, de acordo com o que tiver sido usado nas etapas anteriores 
raster_list = list()

for(year in seq(number_of_projections)){

  filepaths = vector()
  filenames = vector()

  for(i in seq(number_of_instances)){

    myRespName = paste('instancia_', i, sep='')
    filepaths = append(filepaths, paste(projectFolder, '/outputs/projections/', gsub('_', '.', myRespName), '/proj_', year, 'kyr.asc', sep=''))
    filenames = append(filenames, paste('instance_', i, '_proj_', year, 'kyr', sep=''))

  }

  raster_stack = stack(filepaths)
  names(raster_stack) = filenames
  raster_list[[paste('proj', year, sep='')]] =  raster_stack
}

# inspecionando um dos resultados (projecao para 1000 anos atras)
plot(raster_list[['instance_1_proj_1kyr']])

# exemplo de predicao com analise de incerteza para 3 mil anos atras
# CORES: 
#  branco = provavel habitat não adequado; 
#  amerelo = maior incerteza se seria adequado; 
#  laranja = menor incerteza se seria adequado
#  verde escuro = provavel area de habitat adequado

plot( sum(raster_list[['proj3']] > 0.1), col=c('white', 'yellow', 'orange', 'darkgreen'))


## CONTINUAR DAQUI ##
# - rever esquema de amostragem para cobertura dos intervalos de 
#   incerteza dos registros fossies (esta ficando mto grande);
# - rever o objeto com os retultados finais (deixar ealmente uma
#   lista para as n instancias e m projecoes temporais?);
# - ver o q daria pra deixar de salvar, para nao ocupar tanto 
#   memoria fisica.
#####################