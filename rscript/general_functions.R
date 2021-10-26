library(tuneRanger)
library(gstat)
library(reshape2)
library(e1071)
library(caret)
library(randomForest)
library(stringr)
library(stars)
library(sf)
library(dplyr)
library(gdalUtils)
library(raster)
library(automap)
library(fields)
library(interp)
library(mgcv) 
library(purrr)
library(furrr)
library(doParallel)
library(future.apply)
# require(ggforce)
# require(OGC)
# require(RSelenium)
# require(pdftools)
# require(mapview)
# require(leaflet)
# require(RSAGA)
# require(Rsagacmd)
# library(patchwork)
# library(viridis)
# require(ggplot2)
# require(gridExtra)
# require(RefManageR)

which.median <- function(x) which.min(abs(x - median(x)))

getInt <- function() {
  
  conPostgres()
  # import manually categorised data
  inc <- st_read('data/includes.gpkg') %>% 
    # filter(interpolate == T) %>% 
    mutate_if(is.factor,as.character)
  # buffer selected polygons
  inc.buff <- inc %>% st_buffer(250) %>% 
    mutate(a_loc = paste(req,a,sep='_'))
  # write buffered polygons
  st_write(inc.buff,'data/inc_buff.shp',
           delete_dsn=T)
  
  # manually trace / mark out features for interpolation 
  # layer name: 'interpolation'
  
  # import interpolation polygon
  conPostgres()
  int <- st_read(con,c('quarry','interpolate')) %>% 
    left_join(inc %>% st_drop_geometry(), by=c('include_id'='pkey')) %>% 
    mutate(a_locs = paste0(userDataDir,'/',layer,
                           '/diffs/',req,'_',a,
                           '.tif'),
           b_locs = paste0(userDataDir,'/',layer,
                           '/diffs/',req,'_',b,
                           '.tif'))
}

conPostgres <- function() {
  
  # install.packages("RPostgreSQL")
  require("RPostgreSQL")
  
  # create a connection
  # save the password that we can "hide" it as best as we can by collapsing it
  load('private/pw.RData')
  
  # loads the PostgreSQL driver
  drv <- dbDriver("PostgreSQL")
  # creates a connection to the postgres database
  # note that "con" will be used later in each connection to the database
  con <<- dbConnect(drv, dbname = 'mal',
                    host = "localhost", port = 5432,
                    user = 'barneyharris', password = pw)
  rm(pw) # removes the password
}

st_write_withPkey <- function(con,sf,schema='',tname,pkeycol='pkey') {
  # sf <- iw
  st_write(sf,con,c(schema,tname))
  if (schema != '') {
    schema <- paste0(schema,'.')
  }
  dbSendQuery(con,paste0('ALTER TABLE ',schema,tname,
                         ' ADD COLUMN pkey SERIAL PRIMARY KEY;'))
}

clearTempTables <- function(con) {
  lapply(dbListTables(con)[str_detect(dbListTables(con),'_temp$')],
         dbRemoveTable, conn=con)
}

# useful print paste function
pp <- function(...) {
  print(paste0(list(...),collapse=''))
}

'%!in%' <- function(x,y)!('%in%'(x,y))

genSections <- function(con,uidNum,secLength,
                        secDensity=sectionDensity,
                        xDensity=1) {
  print('generating section lines...')
  # original earthwork geometry lines
  ewLinesWgs <- st_read(con, c('rs','earthwork_geometry')) %>% 
    filter(uid == uidNum) %>% st_transform(4326)
  ewLines <- st_read(con, c('rs','earthwork_geometry')) %>% 
    filter(uid == uidNum) %>% st_cast('LINESTRING')
  ewBuff <- st_read(con, c('rs','earthwork_geometry')) %>% 
    filter(uid == uidNum) %>% st_buffer(secLength) %>% 
    st_union()
  # export for raster processing
  ewBuff %>% 
    st_write('/srv/mal_data/rmal/shapefiles/ewbuff_temp.gpkg',
             delete_dsn=T)
  
  leaflet() %>% addTiles() %>% 
    addPolylines(data=st_transform(ewLines,4326))
  # convert lines to points, numbered in sequence
  ewPoints <- st_line_sample(ewLines, density=secDensity,
                             type='regular') %>%
    st_as_sf() %>% st_cast('POINT') %>% mutate(cat = row_number()) %>% 
    arrange_at('cat') %>% st_transform(4326)
  
  # calculate bearings and difference between bearings between points
  ewPointsSp <- ewPoints %>% as_Spatial()
  ewBearings <- bearing(ewPointsSp)
  ewPoints$bearings <- ewBearings
  ewPoints$diff <- diff(c(ewPoints$bearings,NA),lag=1)
  # create directionless bearing for purposes of comparison with objects
  ewPoints$pos_bearing <- (ewPoints$bearings + 180) %% 180
  
  ewPointsOsgb <- st_as_sf(ewPoints) %>% st_transform(27700)
  # generate section lines
  s1 <- as.data.frame(destPoint(ewPointsSp, ewBearings-90, 
                                (secLength/2) ))
  s1Sf <- st_as_sf(s1[complete.cases(s1),],
                   coords=c('lon','lat'), crs=4326) %>% 
    mutate(id = row_number(),
           s = 's1')
  
  st_write(ewPoints, con, 'secpoints_temp')
  print('written section points to db as... secpoints_temp')
  # st_write(s1Sf, con, 'ewpoints_s1_temp')
  s2 <- as.data.frame(destPoint(ewPointsSp, ewBearings+90, 
                                (secLength/2) ))
  s2Sf <- st_as_sf(s2[complete.cases(s2),],
                   coords=c('lon','lat'),
                   crs=4326) %>% 
    mutate(id = row_number(),
           s = 's2')
  
  sectionLines <- rbind(s1Sf,s2Sf) %>% group_by(id) %>% 
    summarise(geometry = st_combine(geometry)) %>%
    st_cast("LINESTRING")
  
  # convert section lines into SP points
  sectionLinesGB <- st_transform(sectionLines, 27700)
  st_write(sectionLinesGB, con, 'seclines_temp')
  print('written section lines to db as... seclines_temp')
  
  sectionP <- st_line_sample(sectionLinesGB, density=xDensity)
  # st_write(sectionP, con, 'sectionpoints_temp')
  sectionPoints <- lapply(sectionP, function(x) {
    co <- st_coordinates(x)
    sp::SpatialPoints(co,CRS(st_crs(27700)$proj4string))
  })
  
  allPoints <- do.call(rbind,sectionPoints)
  
  
  out <- list(lines = sectionLines,
              points = sectionPoints,
              pointsDf = allPoints,
              sf_points = sectionP,
              ew_points = ewPointsOsgb,
              ew_buff = ewBuff)
  return(out)
}

# functions for 1_process_quarry_locations ----

# bufferedPoly <- qGrps$SUSW
# bufferedPoly <- qGrps$STNE
# max request limit = 625455396.147
# x <- 1
# bufferedPoly <- int[x,]
# whichProd <- "LIDAR Point Cloud"
# whichYears <- c(int[x,]$a,int[x,]$b)
# bufferedPoly <- sfFilt
getLidar2 <- function(bufferedPoly,
                      whichProd="LIDAR Tiles DTM",
                      whichYears,
                      minSurvey = 2,
                      userDataDirRoot='tmp',
                      overwrite=T) {
  
  if (is.null(bufferedPoly$tile50k_name_char)) {
    bufferedPoly$tile50k_name_char <- bufferedPoly$layer
  }
  
  # define chrome options
  eCaps <- list(chromeOptions = list(
    args = c(
      '--disable-gpu'
      ,'--headless',
      '--window-size=1280,800'
    )
  ))
  rD <- rsDriver(browser = "chrome",
                 chromever = "81.0.4044.138",
                 extraCapabilities = eCaps,
                 port =
                   as.integer(base::sample(seq(32768,65535, by=1),1)))
  # rD <- RSelenium::rsDriver(
  #   browser = "firefox",
  #   extraCapabilities = list(
  #     "moz:firefoxOptions" = list(
  #       args = list('--headless')
  #     )
  #   ),
  #   port = 
  #     as.integer(base::sample(seq(32768,65535, by=1),1))
  # )
  
  remDr <- rD[["client"]]
  
  browseQuery <- function(remDr,bufferedPoly) {
    bufferedPoly <- st_union(bufferedPoly) %>% st_sf
    st_write(bufferedPoly, dsn=paste0('data/temp.shp'),
             delete_dsn=T)
    simplePoly <- paste0(getwd(),'/data/temp.shp')
    simplePolyZip <- paste0(getwd(),'/data/temp.zip')
    simplePolyFiles <- paste0(getwd(),'/data/temp*')
    # zip shapefiles into archive
    system(paste0('zip ',simplePolyZip,' ',simplePolyFiles))
    
    # navigate to the EA lidar portal
    remDr$navigate("https://environment.data.gov.uk/DefraDataDownload/?Mode=survey")
    # upload zipped shape file
    suppressMessages({
      try(remDr$findElement("id", "fileid"),TRUE)
    })
    while (remDr$status == 7) {
      Sys.sleep(2)
      print('waiting for portal to load...')
      suppressMessages({
        try(remDr$findElement("id", "fileid"),silent=TRUE)
      })
    }
    webElem <- remDr$findElement("id", "fileid")
    webElem$sendKeysToElement(list(simplePolyZip))
    print('uploading shape file...')
    # wait for upload to complete (2 seconds)
    Sys.sleep(5)
    # find 'Get Tiles' button
    getTiles <- remDr$findElement(using = 'css selector', ".grid-item-container")
    # click 'Get Tiles' button
    getTiles$clickElement()
    # sys.sleep ?
    suppressMessages({
      try(remDr$findElement(using = 'css selector', '.data-ready-container'),TRUE)
    })
    
    i <- 0
    while (remDr$status == 7) {
      Sys.sleep(5)
      print(paste0('waiting for tiles to be returned...'))
      suppressMessages({
        try(remDr$findElement(using = 'css selector', '.data-ready-container'),TRUE)
      })
      i <- i + 1
      if (i > 12) {
        print('error with shape file...')
        return('shapefile error')
      }
    }
    print('tiles returned!')
  } # / browseQuery
  browseResult <- browseQuery(remDr,bufferedPoly)
  if (browseResult == 'shapefile error') {
    print(browseResult)
    return(browseResult)
  }
  l <- leaflet() %>% addTiles() %>% 
    addPolygons(data=st_transform(bufferedPoly,4326))
  print(l)
  print('searching available tiles...')
  # select products DTMs
  desiredProds <- whichProd
  # desiredProds <- "LIDAR Tiles DTM"
  # desiredProds <- "LIDAR Point Cloud"
  prodElem <- remDr$findElement(using = 'css selector', '#productSelect')
  prodList <- unique(prodElem$selectTag()$text)
  prodsIndex <- which(prodList %in% desiredProds)
  
  xP <- paste0('//*[@id="productSelect"]/option[',prodsIndex,']')
  webElem <- remDr$findElement(using = 'xpath', 
                               value = xP)
  webElem$clickElement()
  webElem$getElementText()
  # check which year available
  yrElem <- remDr$findElement(using = 'css selector', '#yearSelect')
  yrList <- unique(yrElem$selectTag()$text)
  
  if (desiredProds == "LIDAR Tiles DTM") { 
    # cycle through years, selecting 1m res and recording tiles names
    # x <- 1
    tileList <- lapply(1:length(yrList), function(x) {
      yr <- yrList[x]
      xP <- paste0('//*[@id="yearSelect"]/option[',x,']')
      webElem <- remDr$findElement(using = 'xpath', 
                                   value = xP)
      webElem$clickElement()
      # now cycle through res
      resElem <- remDr$findElement(using = 'css selector', '#resolutionSelect')
      resVec <- unique(resElem$selectTag()$text)
      # pick only 1m
      if (length(which(resVec == 'DTM 1M')) == 0) {
        return(NULL) } else { r <- which(resVec == 'DTM 1M') }
      
      resElem$clickElement() # open drop down
      xP <- paste0('//*[@id="resolutionSelect"]/option[',r,']')
      webElem <- remDr$findElement(using = 'xpath', 
                                   value = xP)
      webElem$clickElement() # select 1m res
      tileLinks <- remDr$findElement(using = 'css selector', '.data-ready-container')
      tileLinks.a <- tileLinks$findChildElements('tag', 'a')
      tiles <- unlist(lapply(tileLinks.a, function(x) x$getElementAttribute('href')))
      
      return(tiles)
    })
    
    # name list by years
    names(tileList) <- yrList
    # remove nulls (years with no 1m res)
    tileList[unlist(lapply(tileList,is.null))] <- NULL
    
  }
  
  if (desiredProds == "LIDAR Point Cloud") {
    # yr <- "2011"
    tileList <- lapply(whichYears, function(yr) {
      x <- which(yrList==yr)
      if (length(x) == 0) {
        print(paste0('year ',yr,' not available as LAZ'))
        return(NULL)
      } 
      xP <- paste0('//*[@id="yearSelect"]/option[',x,']')
      webElem <- remDr$findElement(using = 'xpath', 
                                   value = xP)
      webElem$clickElement()
      
      tileLinks <- remDr$findElement(using = 'css selector', '.data-ready-container')
      tileLinks.a <- tileLinks$findChildElements('tag', 'a')
      tiles <- unlist(lapply(tileLinks.a, function(x) x$getElementAttribute('href')))
      
      return(tiles)
    })
    
    # name list by years
    names(tileList) <- whichYears
    tileList[unlist(lapply(tileList, is.null))] <- NULL
  }
  
  # extract tile names from download URLs
  # x <- names(tileList)[1]
  tileNames <- lapply(names(tileList), function(x) {
    unlist(lapply(str_split(tileList[[x]], 
                            paste0(x,'-')),function(y) substr(y[2],1,6)))
  })
  names(tileNames) <- names(tileList)
  # convert data to tile names with lists of years
  tilesYears <- lapply(unique(unlist(tileNames)), function(tile) {
    allYears <- lapply(names(tileNames), function(year) {
      if (tile %in% tileNames[[year]]) year
    })
    allYears[unlist(lapply(allYears,is.null))] <- NULL
    return(unlist(allYears))
  })
  names(tilesYears) <- unique(unlist(tileNames))
  
  # minimun number of years survey 3, remove the rest
  if (minSurvey > 0) {
    tilesYears[unlist(lapply(tilesYears, function(x) length(x) < minSurvey))] <- NULL
    if (length(tilesYears) == 0) {
      er <- 'no tiles with sequential surveys found...'
      print(er)
      return(er)
    }
  }
  
  allLinks <- as.character(unlist(tileList))
  dlLinks <- allLinks[str_detect(allLinks,paste(names(tilesYears),collapse = '|'))]
  
  # output URLs as list for Wget
  fileName <- paste0(unique(bufferedPoly$tile50k_name_char,'_list.txt'))
  
  write.table(dlLinks,
              file=paste0('wget/',fileName),
              quote = F,row.names=F,col.names = F)
  print(paste0('written download list to ... wget/',fileName))
  
  # close selenium
  remDr$close()
  rD$server$stop()
  gc()
  
  # create folder structure
  folderPath <- paste0(userDataDir,'/',userDataDirRoot,'/',unique(bufferedPoly$tile50k_name_char))
  if (!dir.exists(folderPath)) {
    dir.create(folderPath)
    lapply(unique(unlist(tilesYears)),function(x) dir.create(paste0(folderPath,'/',x)))
  }
  
  # overwrite=T
  if (overwrite) {
    system(paste0('rm ',folderPath,' -R'))
    dir.create(folderPath)
    lapply(unique(unlist(tilesYears)),function(x) dir.create(paste0(folderPath,'/',x)))
  }
  
  # download and uncompress EA lidar with magic!
  # x <- 1
  system(paste0('cat ',getwd(),'/',paste0('wget/',fileName),' | parallel --gnu ',
                shQuote(paste0('wget {} -P ',folderPath))))
  
  # extract to yearly folders 
  yrs <- unique(unlist(tilesYears))
  # x <- yrs[2]
  lapply(yrs, function(x) {
    zips <- paste0(folderPath,'/',
                   list.files(folderPath)[grep(paste0("*(",x,").*zip$"),list.files(folderPath))])
    if (length(zips) > 0) { 
      lapply(zips, function(y) {
        system(paste0("unzip -n ",y,
                      ' -d ',folderPath,'/',x,'/'))
      })
    }
  })
  
  return(folderPath)
}

# function to find a wider area of overlapping DTMs, for use in the UL
# dataset
# folderName <- 'STSE'
moreOverlaps <- function(folderName,
                         folderPath='/home/barneyharris/user_quarry_data/dtms/',
                         minSurvey=2,    # min number of surveys
                         polBuffer=10,   # how big to make dummy polygon
                         patchSize=50,   # size of DEM patch
                         maxDiff=0.1,    # max allowable difference
                         singleThread=T  # F = multicore
) {
  fd <- paste0(folderPath,folderName)
  print(paste0('processing folder...',fd))
  f <- list.files(fd,pattern='*.tif$',
                  full.names = T,
                  recursive = T)
  # remove any tifs from diffs path
  f <- f[!str_detect(f,'diffs')]
  # create index 
  indexLoc <- paste0(fd,'/',basename(fd),'_tindex.shp')
  if (file.exists(indexLoc)) {
    system(paste0('rm ',str_replace(indexLoc,'.shp','*'))) }
  suppressWarnings({gdaltindex(indexLoc,
                               f)})
  tIndex <- st_read(indexLoc,
                    crs=27700,quiet = T,
                    stringsAsFactors = F) %>% 
    mutate(tile_id = 1:nrow(.),
           wkt = st_as_text(geometry)) %>% 
    group_by(wkt) %>% 
    add_tally() %>% 
    dplyr::filter(n >= minSurvey) %>% 
    mutate(year = basename(dirname(location)))
  
  toVrt <- tIndex %>% 
    dplyr::group_by(year) %>% 
    dplyr::summarise(locs = list(location)) %>% 
    mutate(year = as.numeric(year)) %>% 
    dplyr::arrange(year)
  
  # create dir for vrts
  if (!dir.exists(paste0(fd,'/vrts'))) {
    dir.create(paste0(fd,'/vrts'))
  }
  # generate vrts
  # x <- 1
  st <- lapply(1:nrow(toVrt), function(x) {
    vrt <- paste0(fd,'/vrts/',
                  toVrt[x,]$year,'.vrt')
    print(vrt)
    gdalbuildvrt(unlist(toVrt[x,]$locs),
                 output.vrt = vrt
                 # ,tr=c(1,1)
                 ,a_srs=st_crs(27700,parameters=T)$proj4string
    )
    r <- raster(vrt)
    crs(r) <- st_crs(27700,parameters=T)$proj4string
    return(r)
  })
  
  names(st) <- toVrt$year
  
  # now subtract earliest survey from latest survey
  if (!dir.exists(paste0(fd,'/diffs'))) {
    dir.create(paste0(fd,'/diffs'))
  }
  
  co <- c('TILED=YES','COMPRESS=DEFLATE')
  print('generating difference rasters...')
  
  if (singleThread == F) {
    w <- max(1:(length(toVrt$year)-1))
    if (w > 6) w <- 6
    plan(multisession, workers=w)
  } else w <- 1
  # w <- 1
  # x <- 2
  diffs <- 1:(length(toVrt$year)-1) %>% 
    future_map(.f = function(x) {
      minYear <- as.character(toVrt$year[x])
      maxYear <- as.character(toVrt$year[x+1])
      print(paste0('processing ',minYear,' ',maxYear))
      diffRas <- tryCatch({st[[maxYear]] - st[[minYear]]},
                          error = function(x) { return('error in rasters..different res?') })
      # check raster returned
      if (is.Raster(diffRas)) {
        # check it's not all NA
        if (freq(diffRas, value=NA)/length(diffRas)<1) {
          # check enough non-NA cells for patch size
          if (cellStats(!is.na(diffRas),sum) > (patchSize*4)^2) {
        
            
        # diff
        # diffOut <- paste0(fd,'/diffs/',
        #                   maxYear,'_',minYear,'_diff.tif')
        # writeRaster(diffRas, diffOut, options=co, overwrite=T)
        # diffBinaryOut <- paste0(fd,'/diffs/',
        #                   maxYear,'_',minYear,'_diff_binary.tif')
        # writeRaster(diffRas > -Inf, diffBinaryOut, options=co, overwrite=T)
        # years
        maxOut <- paste0(fd,'/diffs/',maxYear,'.tif')
        minOut <- paste0(fd,'/diffs/',minYear,'.tif')

        # gdaladdo(diffOut, levels=c('8','16','32','64','128'))
        # gdaladdo(maxOut, levels=c('8','16','32','64','128'))
        # gdaladdo(minOut, levels=c('8','16','32','64','128'))
        # make patch grid
        
        diffStars <- st_as_stars(diffRas)
        diffPol <- terra::as.polygons(as(diffRas > -Inf, "SpatRaster")) %>% 
          st_as_sf
        
        diffPatch.sf <- diffPol %>% 
          st_make_grid(.,cellsize=patchSize) %>% 
          st_sf %>% 
          st_join(.,st_buffer(diffPol,-patchSize), st_covered_by) %>% 
          tidyr::drop_na() %>% 
          st_sf
        
        # plot(diffPol)
        # plot(diffPatch.sf$geometry,add=T)
        
        diffPatch.sf$diffThresh <- 1:nrow(diffPatch.sf) %>% 
          map_chr(.f = function(x) {
            isDiff <- any(abs(diffStars[diffPatch.sf[x,]][,,][[1]]) > maxDiff,
                          na.rm=T)
          })
        
        # plot(diffPatch.sf[,'diffThresh'],
        #      main=paste0(maxYear,' - ',minYear,
        #                  ' with diff of ',maxDiff,
        #                  'm or less'))
        
        diffPatch.filt <- diffPatch.sf %>% dplyr::filter(diffThresh == 'FALSE')
        
        if (nrow(diffPatch.filt) == 0) return(NULL)
        
        diffPols <- diffPatch.filt$geometry %>% map(.f = function(g) {
          st_sample(g,3) %>% 
            st_sf %>% 
            mutate(id = 1) %>% 
            group_by(id) %>% st_union() %>% st_cast("LINESTRING") %>% 
            smoothr::smooth('chaikin') %>% st_buffer(polBuffer/2) %>% st_sf
        }) %>% do.call(rbind, .) %>% 
          dplyr::mutate(fid = 1:nrow(.),
                        a_loc = maxOut,
                        b_loc = minOut,
                        a = str_replace(basename(maxOut),'.tif',''),
                        b = str_replace(basename(minOut),'.tif',''))
        
        return(diffPols)
          }
        }
      }
    }) %>% compact()
  
  # write neccessary tifs
  y <- diffs %>% map(~unique(c(.x$a,.x$b))) %>% 
    unlist() %>% unique(.)
  print('writing tifs...')
  y %>% walk(~writeRaster(st[[.x]], paste0(fd,'/diffs/',.x,'.tif'), 
                          options=co, overwrite=T))
  
  return(diffs)
}

# bufferedPoly <- qGrps$STSW
# folderPath <- fp

# bufferedPoly <- int[x,]
# whichProd <- "LIDAR Point Cloud"
# whichYears <- c(int[x,]$a,int[x,]$b)
processLidar <- function(polyGrps=qGrps,folderName,removeDls=T) {
  print(paste0('processing folder...',folderName))
  folderPath <- paste0(userDataDir,'/',folderName)
  bufferedPoly <- polyGrps[[basename(folderPath)]]
  # check which raster tiles overlap with one another...
  # list all files
  f <- list.files(folderPath,pattern='*.tif$',
                  full.names = T,
                  recursive = T)
  # remove any tifs from diffs path
  f <- f[!str_detect(f,'diffs')]
  # create index 
  indexLoc <- paste0(folderPath,'/',basename(folderPath),'_tindex.shp')
  if (file.exists(indexLoc)) {
    system(paste0('rm ',str_replace(indexLoc,'.shp','*'))) }
  suppressWarnings({gdaltindex(indexLoc,
                               f)})
  # build tile index of rasters, group by geometry and
  # discard any tiles which do not overlap more than
  # minimum survey value
  reqSng <- bufferedPoly %>% st_cast('POLYGON') %>% 
    mutate(req_id = 1:nrow(.))
  
  tIndex <- st_read(indexLoc,
                    crs=27700,quiet = T,
                    stringsAsFactors = F) %>% 
    mutate(tile_id = 1:nrow(.),
           wkt = st_as_text(geometry)) %>% 
    group_by(wkt) %>% 
    add_tally() %>% 
    dplyr::filter(n >= minSurvey) %>% 
    mutate(year = basename(dirname(location))) %>% 
    st_join(reqSng, left = F) # of the remaining polygons, keep 
  # only those which intersect with with quarries request polygon
  
  if (nrow(tIndex) == 0) return('no rasters intersecting with quarries...')
  
  # generate VRT lists, grouping rasters first by request polygon ID,
  # then year of survey
  toVrt <- tIndex %>% group_by(req_id,year) %>% 
    summarise(locs = list(location)) %>% 
    mutate(year = as.numeric(year))
  # create dir for vrts
  if (!dir.exists(paste0(folderPath,'/vrts'))) {
    dir.create(paste0(folderPath,'/vrts'))
  }
  # generate vrts
  x <- 1
  st <- lapply(1:nrow(toVrt), function(x) {
    vrt <- paste0(folderPath,'/vrts/',
                  'req',toVrt[x,]$req_id,'_',
                  toVrt[x,]$year,'.vrt')
    print(vrt)
    gdalbuildvrt(unlist(toVrt[x,]$locs),
                 output.vrt = vrt
                 # ,tr=c(1,1)
                 ,a_srs=st_crs(bufferedPoly,parameters=T)$proj4string
    )
    r <- raster(vrt)
    crs(r) <- st_crs(bufferedPoly,parameters=T)$proj4string
    return(r)
  })
  
  names(st) <- paste0('req',toVrt$req_id,'_',
                      toVrt$year)
  
  # now subtract earliest survey from latest survey
  if (!dir.exists(paste0(folderPath,'/diffs'))) {
    dir.create(paste0(folderPath,'/diffs'))
  }
  
  co <- c('TILED=YES','COMPRESS=DEFLATE')
  print('generating difference rasters...')
  x <- 12
  diffsResult <- lapply(unique(toVrt$req_id), function(x) {
    print(paste0('processing request poly... ',x))
    df <- toVrt %>% as.data.frame() %>% 
      filter(req_id == x)
    maxYear <- paste0('req',x,'_',max(df$year))
    minYear <- paste0('req',x,'_',min(df$year))
    diffRas <- tryCatch({st[[maxYear]] - st[[minYear]]},
                        error = function(x) { return('error in rasters..different res?') })
    if (!is.Raster(diffRas)) return(diffRas)
    
    if (!is.na(raster::maxValue(diffRas)) && raster::maxValue(diffRas) > 2) {
      diffOut <- paste0(folderPath,'/diffs/',
                        'req',x,'_',max(df$year),
                        '_',min(df$year),'_diff.tif')
      writeRaster(diffRas, diffOut, options=co, overwrite=T)
      maxOut <- paste0(folderPath,'/diffs/',
                       'req',x,'_',max(df$year),'.tif')
      minOut <- paste0(folderPath,'/diffs/',
                       'req',x,'_',min(df$year),'.tif')
      writeRaster(st[[maxYear]], maxOut, options=co, overwrite=T)
      writeRaster(st[[minYear]], minOut, options=co, overwrite=T)
      gdaladdo(diffOut, levels=c('8','16','32','64','128'))
      gdaladdo(maxOut, levels=c('8','16','32','64','128'))
      gdaladdo(minOut, levels=c('8','16','32','64','128'))
      return('succesfully produced differences rasters')
    } else return('no differences over 2 metres')
  })
  
  # remove original
  if (removeDls) {
    system(paste0('rm ',folderPath,'/*.zip'))
    yearDirs <- setdiff(list.dirs(folderPath),folderPath)
    system(paste0('rm ',paste(yearDirs[!str_detect(yearDirs,'diffs')],
                              collapse=' '),' -R'))
  }
  # plot(r)
  return(diffsResult)
}


# sfPoly <- int[1,]
processLaz <- function(sfPoly) {
  require(lidR)
  print(paste0('processing folder...',sfPoly$layer))
  folderPath <- paste0(userDataDir,'/laz/',sfPoly$layer)
  lazDirs <- list.dirs(list.dirs(folderPath,
                                 recursive = F,
                                 full.names =T),
                       recursive = F,
                       full.names = T)
  
  lapply(lazDirs, function(x) {
    lazPath <- paste0(x,'/')
    ctg <- readLAScatalog(lazPath)
    las_check(ctg)
    sfBuff <- sfPoly %>% st_buffer(50)
    sfDiff <- st_difference(sfBuff,sfPoly)
    
    processPoly <- function(sf, tag) {
      st_write(sf,paste0('data/lidar_poly_',tag,'.shp'),
               delete_dsn=T)
      ctgRoi <- shapefile(paste0('data/lidar_poly_',tag,'.shp'))
      roi <- clip_roi(ctg, ctgRoi)
      
      ras <- grid_terrain(roi, algorithm = tin())
      crs(ras) <- st_crs(27700,parameters=T)$proj4string
      outRoot <- str_replace_all(str_replace(lazPath,paste0(userDataDir,'/laz/'),''),
                                 '/','_')
      outName <- paste0(userDataDir,'/int/',outRoot,sfPoly$req,
                        '_',tag)
      writeRaster(ras,filename = paste0(outName,'.tif'),
                  overwrite=T)
      rasHs <- hillShade(terrain(ras,opt="slope"),
                         terrain(rasDiff,opt="aspect"),
                         filename=paste0(outName,'_hs.tif'),
                         overwrite=T)
    }
    
    try(processPoly(sfBuff,'buff'))
    try(processPoly(sfDiff,'diff'))
    
  })
  
}


processDiffs <- function(fname,minDiff=0.5,
                         minArea=500,
                         maxArea=10000,
                         refresh=F) {
  pp('processing folder...',fname)
  folderPath <- paste0(userDataDir,'/',fname)
  diffPath <- paste0(folderPath,'/diffs')
  # load rasters
  diffs <- lapply(list.files(diffPath,full.names = T,pattern='*diff.tif'),raster)
  if (length(diffs) == 0) return('no difference rasters in folder..')
  # calc function to remove cells less than n difference
  pp('processing difference rasters with minimum difference of... ',minDiff)
  # rCalc <- function(x) { x[dplyr::between(x,-minDiff,minDiff)] <- NA; return(x) }
  rCalc <- function(x) { x[x < minDiff] <- NA; return(x) }
  diffsCalc <- lapply(diffs, calc, rCalc)
  
  # sieve rasters, removing clumps smaller than n
  # x <- diffsCalc[[1]]
  pp('sieving rasters to leave only areas larger than ',minArea,' sq. m')
  pp('sieving rasters to leave only areas smaller than ',maxArea,' sq. m')
  diffsSieve <- lapply(diffsCalc, function(x) {
    clump <- raster::clump(x)
    clumpDf <- as.data.frame(freq(clump))
    # put these into a vector of clump ID's to be removed
    excludeID <- clumpDf$value[which(clumpDf$count < minArea)]
    excludeID <- c(excludeID,clumpDf$value[which(clumpDf$count > maxArea)])
    clump.sieve <- clump
    clump.sieve[clump %in% excludeID] <- NA
    # plot(clump.sieve)
    return(clump.sieve)
  })
  
  # vectorize and add to gpkg
  pp('adding polygons to vector layer...')
  pols <- do.call(rbind,lapply(1:length(diffsSieve), function(x) {
    pol <- st_simplify(st_as_sf(st_as_stars(diffsSieve[[x]], crs=st_crs(27700)),
                                as_points = F, merge=T,10,
                                crs=27700)) %>% 
      mutate(ras = basename(diffs[[x]]@file@name) )
    st_crs(pol) <- st_crs(27700)
    return(pol)
  }))
  
  if (nrow(pols) > 0) {
    st_write(pols,paste0('data/aoi_poly_wfiles.gpkg'),
             layer=fname,
             append = F)
  }
  
  # require(ggplot2)
  # a <- as.data.frame(diffs[[1]],xy=T)
  # b <- as.data.frame(diffsSieve[[1]],xy=T)
  # plot(diffsSieve[[1]])
  # ggplot() +
  #   # geom_raster(data = a, aes(x = x, y = y, fill=layer)) + 
  #   geom_raster(data = b, aes(x = x, y = y, fill=clumps)) + 
  #   coord_quickmap()
  
}

# functions for 2_interpolation and analyses  ----

# this is a modified version of GSIF:buffer.dist, used in the RF interpolation
# method. It has as additional check for finding classes without points.

# observations <- trainingSp["elev"]
# predictionDomain <- gridPix[1]
# classes=classesElev # groups elev values
# width=maxSize
buffer.dist2 <- function(observations, predictionDomain, classes, width, ...) {
  if(missing(width)){ width <- sqrt(areaSpatialGrid(predictionDomain)) }
  if(!length(classes)==length(observations)){ stop("Length of 'observations' and 'classes' does not match.") }
  ## remove classes without any points:
  xg = summary(classes, maxsum=length(levels(classes)))
  selg.levs = attr(xg, "names")[xg > 0]
  
  if(length(selg.levs)<length(levels(classes))){
    fclasses <- as.factor(classes)
    fclasses[which(!fclasses %in% selg.levs)] <- NA
    classes <- droplevels(fclasses)
  }
  
  ## derive buffer distances
  s <- list(NULL)
  
  for(i in 1:length(levels(classes))) {
    if (nrow(observations[which(classes==levels(classes)[i]),1]) > 0) {
      
      s[[i]] <- raster::distance(
        rasterize(
          observations[which(classes==levels(classes)[i]),1]@coords, y=raster(predictionDomain)),
        width=width)
    } else print(i)
  }
  s <- s[sapply(s, function(x){!is.null(x)})]
  s <- brick(s)
  s <- as(s, "SpatialPixelsDataFrame")
  s <- s[predictionDomain@grid.index,]
  return(s)
}

# buffFid <- 5
generateCoVars <- function(buffFid) {
  
  demFolder <- file.path(userDataDir,'edina','terrain-5-dtm_3791534')
  intBuff <- st_read('data/interpolation_polygons_500m_buff.kml') %>% 
    st_transform(27700)
  
  gdaltindex(index_file = paste0(demFolder,'/tindex.shp'),
             gdal_file = list.files(demFolder,pattern='*.asc$',
                                    full.names = T,
                                    recursive = T))
  tIndex <- st_read(paste0(demFolder,'/tindex.shp'),
                    crs=27700,quiet = T,
                    stringsAsFactors = F)
  # x <- 5
  callGrass <- function(x) {
    ras <- tIndex %>% st_join(intBuff[intBuff$fid==x,],left=F) %>% 
      dplyr::select(location)
    gdalbuildvrt(ras$location,'raster/temp_vrt.vrt',overwrite = T)
    r <- raster('raster/temp_vrt.vrt')
    
    rAg <- aggregate(r, 4)
    writeRaster(rAg, 'raster/temp_ag.tif',overwrite=T)
    # grass python script cam be editted in RStudio!
    pyLoc <- paste0(getwd(),'/python/processDEM.py')
    
    intID <- paste0(x,'_fid')
    
    # geomorphon search radius (in metres)
    geosearch <- '500'
    # run grass via system command, executing above python script
    system2('grass',
            paste(shQuote(gdb),
                  shQuote('--exec'),
                  shQuote(pyLoc),
                  shQuote(file.path(getwd(),'raster/temp_ag.tif')),
                  shQuote(file.path(userDataDir,'coarse')),
                  shQuote(intID),
                  shQuote(geosearch)
            )
            # ,stderr = paste0(getwd(),'/logs/pyfunc_err.txt') # comment out to see output
            ,stdout = paste0(getwd(),'/logs/processDEM_stdout.txt') # comment out to see output
    )
    
    grassOut <- read.table('logs/processDEM_stdout.txt',sep='\n',
                           stringsAsFactors = F) %>% slice_tail()
    grassOut <- str_split(grassOut$V1,', ')[[1]]
    grassFiles <- paste0(userDataDir,'/','coarse/',grassOut,'_',intID,'.tif')
    
    grassRas <- lapply(grassFiles, raster)
    
    grassStars <- read_stars(unlist(grassFiles))
    return(list(ras = grassRas,
                st = grassStars))
    
  } # per interpolation polygon
  
  cov <- callGrass(buffFid)
  return(cov)
}

# dummy data generator
dummyData <- function() {
  s <- 10
  xy <- do.call(rbind,lapply(-s:s, function(x) { data.frame(x = x, y=-s:s) } ))
  xy$z <- sqrt(abs((xy$x*xy$x)+(xy$y*xy$y)))
  
  # df <- xy2
  foursUp <- function(df) {
    df1 <- df
    df2 <- df
    df3 <- df
    # plus x and y the same
    df1[,'x'] <- df1[,'x'] + s*2
    # plus x and minus y
    df2[,'x'] <- df2[,'x'] + s*2
    df2[,'y'] <- df2[,'y'] - s*2
    # x the same and minus y
    df3[,'y'] <- df3[,'y'] - s*2
    bind_rows(df,df1,df2,df3)
  }
  
  xy2 <- foursUp(xy)
  xy3 <- foursUp(xy2)
  xy4 <- foursUp(xy3)
  xyRas <- rasterFromXYZ(xy4)
  plot(xyRas)
  writeRaster(xyRas,'raster/dummy.tif',overwrite=T)
  xyStars <- read_stars('raster/dummy.tif')
  smpMask <- st_as_sf(xy4[sample(1:nrow(xy4),3),],
                      coords = c('x','y')) %>% 
    mutate(id = 1) %>% 
    group_by(id) %>% st_union() %>% st_cast("LINESTRING") %>% 
    smoothr::smooth('chaikin') %>% st_buffer(s/2) %>% st_sf
  
  a_ras <- xyStars
  b_ras <- xyStars
  pol <- smpMask
  st_crs(pol) <- st_crs(27700)
  st_crs(a_ras) <- st_crs(27700)
  st_crs(b_ras) <-st_crs(27700)
  plot(a_ras)
  plot(pol, add=T)
}

# hillshade function
st_hillshade <- function(st_ras) {
  slope.aspect <- terrain(as(st_ras, "Raster"),opt=c('slope','aspect'))
  hs <- hillShade(slope.aspect$slope,slope.aspect$aspect)
  writeRaster(hs,'raster/hs.tif',overwrite=T)
  hsSt <- read_stars('raster/hs.tif')
  return(hsSt)
}

# st_ras <- a_ras
st_slopeaspect <- function(st_ras, formatOut='Stars') {
  st_r <- as(st_ras, "Raster")
  slope.aspect <- terrain(st_r,opt=c('slope','aspect'))
  st_r <- mask(st_r,slope.aspect$slope)
  writeRaster(st_r,'raster/elev.tif',overwrite=T)
  
  writeRaster(slope.aspect$slope,'raster/slope.tif',overwrite=T)
  writeRaster(slope.aspect$aspect,'raster/aspect.tif',overwrite=T)
  
  if (formatOut == 'Stars') {
    output <- read_stars(c('raster/elev.tif','raster/slope.tif','raster/aspect.tif'))
  } else output <- brick(st_r,slope.aspect$slope,slope.aspect$aspect)
  return(output)
}
# pol <- int[1,]
# pol <- polygonDf[1,]
# offsetPol <- F
# bufferDem <- 10
loadLidar <- function(pol, 
                      bufferDem=50, # how much DEM around hole to include?
                      offsetPol=F,
                      offsetDir=1,
                      plotLidar=T) {
  # buffer around interpolation polygon to select ras LIDAR tiles
  b <- bufferDem
  
  # offset by buffer amount
  if (offsetPol) {
    pol_bbox <- st_geometry(pol) %>% st_bbox()
    pol_width <- pol_bbox$xmax-pol_bbox$xmin
    shiftValue <- b - pol_width
    if (offsetDir == 1) offdir <- c(+shiftValue,+shiftValue)
    if (offsetDir == 2) offdir <- c(+shiftValue,-shiftValue)
    if (offsetDir == 3) offdir <- c(-shiftValue,+shiftValue)
    if (offsetDir == 4) offdir <- c(-shiftValue,-shiftValue)
    
    pol.mod <- pol
    g <- st_geometry(pol.mod) + 
      matrix(data = offdir, ncol = 2)
    st_geometry(pol.mod) <- g
    st_crs(pol.mod) <- st_crs(27700)
    # plot(pol$geom %>% st_buffer(500))
    # plot(pol,add=T)
    # plot(pol.mod$geom,add=T)
    pol <- pol.mod
  }
  
  # load rasters into memory, and clip by polygon
  # refCol <- 'a'
  loadRas <- function(pol,b,refCol) {
    # for us with ref col 'a'
    tIndex <- st_read(paste0(userDataDir,'/dtms/',pol$layer,'/',
                             pol %>% st_drop_geometry %>%
                               dplyr::select(all_of(refCol))
                             ),
                      quiet=T
                      )

    tiles <- pol %>% st_join(tIndex) %>% mutate_if(is.factor,as.character)
    locsExist <- rasFiles <- tiles$location[
      tiles$location %>%  map(~file.exists(.x)) %>% unlist()
      ]
    if (length(locsExist) != length(tiles$location)) {
      print('warning...missing tiles detected!') }
    if (length(locsExist) > 1) {
      t <- lapply(rasFiles, raster)
      t$fun <- mean
      r <- do.call(mosaic, t)
      rStars <- st_as_stars(r)
    } else r <- raster(locsExist)
    
    polBuff <- pol %>% st_buffer(b)
    polBuff.r <- fasterize::fasterize(polBuff, raster(r))
    
    # checks whether raster extents overlap
    overlap <- !is.nan(cellStats(r - polBuff.r,mean))
    if (!overlap) return(NULL)
    
    clippedRas.r <- trim(mask(r, polBuff.r))
    writeRaster(clippedRas.r,paste0('raster/',refCol,'.tif'),
                overwrite=T)
    clippedRas <- read_stars(paste0('raster/',refCol,'.tif'))
    
    return(clippedRas)
  }
  
  # load LIDAR rasters as stars objects
  
  a_ras <- loadRas(pol, b, 'a') # latest survey
  b_ras <- loadRas(pol, b, 'b') # earliest survey
  
  if (plotLidar) {
    
    gridExtra::grid.arrange(
      ggplot() + 
        geom_stars(data=a_ras) + 
        geom_sf(data=pol,fill=NA,color='white') + 
        coord_sf(datum = sf::st_crs(27700)) + 
        theme_bw() + 
        theme(legend.position = 'none',
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) + 
        ggtitle(paste0('DEM A: ',pol$a)),
      ggplot() + 
        geom_stars(data=b_ras) + 
        geom_sf(data=pol,fill=NA,color='white') + 
        coord_sf(datum = sf::st_crs(27700)) + 
        theme_bw() + 
        theme(legend.position = 'none',
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) + 
        ggtitle(paste0('DEM B: ',pol$b)),
      nrow=1
    )
  }
  
  return(list(a = a_ras,
              b = b_ras,
              pol = pol))
}


# ras <- tiles$a
# cutpoly <- tiles$pol
# sampleRasPer <- 0
processTiles <- function(ras, 
                         cutpoly=NULL, # cut hole with polygon
                         increaseHole=0, # adjust size of hole from polygon within DEM
                         sampleRasPer=0, # randomly select percent of raster e.g. 0.1
                         addOGC=T,
                         addSlopeAspect=T) {
  
  # # generate slope aspect as raster brick
  # ras.sa.r <- st_slopeaspect(ras,formatOut ='Raster')
  # #
  # # # generate oblique coordinates
  # ras.ogc.r <- makeOGC(as(ras, "Raster"), 8)
  # ras.r <- stack(ras.sa.r,ras.ogc.r)
  
  ras.r <- as(ras, "Raster")
  names(ras.r)[1] <- 'elev'
  
  # cookie cut interplation mask polygon A raster and convert
  # resulting raster to SF points for interpolation
  if (!is.null(cutpoly)) {
    if (increaseHole > 0) cutpoly <- cutpoly %>% st_buffer(increaseHole)
    
    # convert hole polygon into raster
    pol.r <- fasterize::fasterize(cutpoly, raster = raster(ras.r))
    # cookie cut
    train.r <- mask(ras.r,pol.r,inverse=T)
    test.r <- mask(ras.r,pol.r)
    
    # generate distance stat for hole cells, add to stack
    # ras.r <- stack(ras.r,mask(raster::distance(ras.r$elev),pol.r))
    # names(ras.r)[which(names(ras.r)=='layer')] <- 'hole_distance'
  }
  
  if (sampleRasPer > 0) {
    # stars implementation - slow
    ras.sp <- st_as_sf(ras,as_points = T,merge=F) %>% 
      as_Spatial()
    ras.sp.s <- 
      ras.sp[sample(1:nrow(ras.sp),
                    ceiling(nrow(ras.sp)*sampleRasPer),
                    replace = F),]
    
    ras.r.s <- dropLayer(rasterize(ras.sp.s, as(ras, 'Raster')),'ID')
    # incorporate cut poly hole, if produced
    if (!is.null(cutpoly)) ras.r.s <- mask(ras.r.s,train.r)
    train.r <- mask(ras.r,ras.r.s)
    test.r <- mask(ras.r,ras.r.s,inverse=T)
  }
  
  out <- lapply(list(train = train.r,
                     test = test.r,
                     all = ras.r),
                function(x) {
                  ras.st <- st_as_stars(x)
                  st_crs(ras.st) <- st_crs(27700)
                  return(ras.st)
                })
  
  convertStars <- function(r) {
    r.sf <-  st_as_sf(r, merge=F, as_points = T) %>%
      dplyr::rename(elev = 1) %>% na.omit() # simple feature points
    
    # r.sf <- 
    #   st_as_sf(raster::rasterToPoints(as(r, 'Raster'),spatial=T)) %>% 
    #   dplyr::rename(elev = 1) %>% na.omit() # simple feature points
    
    st_crs(r.sf) <- st_crs(27700)
    r.df <-  st_drop_geometry(r.sf) %>% 
      cbind(st_coordinates(r.sf)) # data.frame
    
    r.sp <-  as(r.sf, "Spatial") # sp points
    r.sp@proj4string <- CRS('+init=epsg:27700')
    r.ras = as(r, "Raster") # raster
    crs(r.ras) <- CRS('+init=epsg:27700')
    return(list(sf = r.sf,
                df = r.df,
                sp = r.sp,
                ras = r.ras))
  }
  # convert data
  
  
  out$train <- convertStars(out$train)
  out$all <- convertStars(out$all)
  
  return(out)
}

# mode = 1: process prepared raster as is and compare interpolated
# surface to both B raster and A raster (itself)

# trainingRas <- foldA$train
# testRas <- foldA$all
# cv <- T
# 
# pd <- prepData[[26]]
# paramData <- cvGrids
# trainingData <- pd$foldA$train
# testData <- pd$foldA$all
# maskPoly = pd$pol
# testCV <- F
# gMap <- 'test'
# gLoc='/home/barneyharris/user_quarry_data/grassdb/quarry/'
# intMethods = 'gspline'
# outputDir = '/media/mal/working_files/quarry/'
# intMethods <- 'nn'

interpolateRas <- function(pd,
                           maskPoly, 
                           paramData,
                           testCV = T,
                           outputDir = '/media/mal/working_files/quarry/',
                           outputTag='testnew',
                           gLoc='/home/barneyharris/user_quarry_data/grassdb/quarry/',
                           # gMap='q1',
                           intMethods=c('nn','idw','ok',
                                        'rfsp','tin','gspline',
                                        'gbicubic','gfilter') ) {
  trainingData <- pd$foldA$train
  testData <- pd$foldA$all
  # object for storing run times
  intTimes <- list()
  # object for storing interpolated rasters
  rasterlist <- list()
  # interpolation models ----
  
  # no parameters needed to be trialled for RF and TIN
  # random forest model, uses buffer distances----
  fit_RF <- function(trainingSp, testSp, tag='',subSample=NULL,
                     classInt=0.01) {
    # save crs from training data
    spCrs <- crs(trainingSp)
    
    # sample hole points
    if (!is.null(subSample)) {
      trainingSp <- trainingSp[sample(nrow(trainingSp),subSample),]
    }
    
    intTemplateCopy <- testSp
    # convert full a_ras to grid of locations to predict elev distances against
    gridPix <- as(testSp, "SpatialPixelsDataFrame")
    
    # split elevation into classes to reduce comp time
    classesElev <- cut(trainingSp$elev,
                       breaks=seq(min(trainingSp$elev),
                                  max(trainingSp$elev),
                                  by=classInt))
    print(paste0('splitting elevation values by ',
                 classInt, ' into ',length(levels(classesElev)),
                 ' classes...'))
    print('calculating buffer distances...')
    maxSize <- max(ncol(testData$ras),nrow(testData$ras))^2
    
    # grid.dist0 <- GSIF::buffer.dist(trainingSp["elev"], 
    #                                 gridPix[1], 
    #                                 # classes=as.factor(1:length(unique(trainingSp$elev))), # all elevs
    #                                 classes=classesElev,# groups elev values
    #                                 width=maxSize # max search radius
    # )
    # 
    grid.dist0 <- buffer.dist2(trainingSp["elev"],
                                    gridPix[1],
                                    # classes=as.factor(1:length(unique(trainingSp$elev))), # all elevs
                                    classes=classesElev,# groups elev values
                                    width=maxSize # max search radius
    )
    
    buffDists <- cbind(grid.dist0@data,gridPix@data)
    grid.dist0@data <- buffDists
    
    ov.elev <- over(trainingSp["elev"], grid.dist0)
    
    
    # find best mtry value
    print('tuning Random Forest to find best "mtry" value...')
    x <- ov.elev[,setdiff(names(ov.elev),'elev')]
    y <- ov.elev[,'elev']
    dn0 <- paste(setdiff(names(ov.elev),c('elev','slope','aspect')),
                 collapse="+")
    fm0 <- as.formula(paste("elev ~ ", dn0))
    
    # tune and construct RF
    m.elev <- tuneRanger::tuneMtryFast(fm0, 
                                       data = ov.elev,
                                       num.treesTry = 250,
                                       stepFactor = 1.5, 
                                       improve = 1e-5,
                                       trace = TRUE, 
                                       plot = FALSE,
                                       doBest = T,
                                       quantreg=TRUE)
    
    elev.rfd <- predict(m.elev, grid.dist0@data, type="quantiles")$predictions
    
    intTemplateCopy$pred_elev <- elev.rfd[,2]
    
    outList <- list(rfSp = rasterFromXYZ(cbind(intTemplateCopy@coords,
                                               intTemplateCopy$pred_elev),
                                         crs=spCrs))
    
    return(outList)
  }
  
  if ('rfsp' %in% intMethods) {
    st <- Sys.time()
    interp_RF <- fit_RF(trainingSp = trainingData$sp,
                        testSp = testData$sp,
                        tag='classed')
    
    intTimes$RF <- Sys.time()-st
    rasterlist$`Random Forest SP` <- list(interp_RF$rfSp)
  }
  
  if ('tin' %in% intMethods) {
    # Triangular Irregular Surface
    st <- Sys.time()
    fit_TIN <- interp::interp( # using {interp}
      x = trainingData$df$X,     # the function actually accepts coordinate vectors
      y = trainingData$df$Y,
      z = trainingData$df$elev,
      xo = testData$df$X, # here we already define the target grid
      yo = testData$df$Y,
      output = "points") %>%
      bind_cols()

    interp_TIN <- raster::rasterFromXYZ(fit_TIN,
                                        crs = crs(trainingData$sf))
    intTimes$TIN <- Sys.time()-st

    rasterlist$`Triangular Irregular Surface` <- list(interp_TIN)
  }
  
  # parameterised interpolation models
  # process CV grids to remove bad combinations (nmin vals > nmax vals)
  # and add run number, add polygon id,
  # test routine ensures parameter with greatest impact on computation time is
  # evenly sampled.
  # df <- paramData$nn
  # paramData <- cvGrids
  
  paramData.c <- 
    paramData %>% 
    map2(.y = names(paramData), 
         .f = function(df, y) {
           if (any(str_detect(names(df),'nmaxVals')) && 
               any(str_detect(names(df),'nminVals'))) {
             df <- df %>% filter(nminVals < nmaxVals) }
           df <- df %>% 
             mutate(run_no = 1:nrow(.),
                    intpol_fid = maskPoly$fid,
                    int_method = y)
           if (testCV) {
             if (y == 'nn' | y == 'idw' | y == 'ok') {
               df <- rbind(df[which.max(df$nmaxVals),],
                           df[which.min(df$nmaxVals),],
                           df[which.median(df$nmaxVals),]
               ) }
             if (y == 'gspline') {
               df <- rbind(df[which.max(df$nminVals),],
                           df[which.min(df$nminVals),],
                           df[which.median(df$nminVals),]
               ) }
             if (y == 'gbicubic') {
               df <- rbind(df[which.max(df$stepVals),],
                           df[which.min(df$stepVals),],
                           df[which.median(df$stepVals),]
               ) }
             if (y == 'gfilter') {
               df <- rbind(df[which.max(df$radVals),],
                           df[which.min(df$radVal),],
                           df[which.median(df$radVals),]
               ) }
           }
           return(df)
         })
  
  if ('nn' %in% intMethods) {
    interp_NNs <- list()
    for (y in 1:length(paramData.c$nn$run_no)) {
      x <- paramData.c$nn$run_no[y]
      pdata <- paramData.c$nn %>% 
        filter(run_no == x)
      print(pdata)
      # Nearest neighbour
      st <- Sys.time()
      fit_NN <- gstat::gstat( # using package {gstat} 
        formula = elev ~ 1,    
        data = trainingData$sp,
        nmax = pdata$nmaxVals,
        nmin = pdata$nminVals,# Number of neighboring observations used for the fit
      )
      interp_NNs[[y]] <- raster::interpolate(testData$ras[[1]], fit_NN)
      tdiff <- Sys.time()-st
      t <- list(val = tdiff,
                unit_chr = units(tdiff))
      intTimes$NN[[y]] <- t
    }
    rasterlist$`Nearest Neighbor` <- interp_NNs
  }
  
  if ('idw' %in% intMethods) { 
    interp_IDWs <- list()
    for (y in 1:length(paramData.c$idw$run_no)) {
      x <- paramData.c$idw$run_no[y]
      pdata <- paramData.c$idw %>% 
        filter(run_no == x)
      
      # Inverse Distance Weighting
      st <- Sys.time()
      fit_IDW <- gstat::gstat( # The setup here is quite similar to NN
        formula = elev ~ 1,
        data = trainingData$sp,
        nmax = pdata$nmaxVals,
        nmin = pdata$nminVals,
        set = list(idp = pdata$idpVals)) # inverse distance power
      
      interp_IDWs[[y]] <- raster::interpolate(testData$ras[[1]], fit_IDW)
      tdiff <- Sys.time()-st
      t <- list(val = tdiff,
                unit_chr = units(tdiff))
      intTimes$IDW[[y]] <- t
    }
    rasterlist$`Inverse Distance Weighted` <- interp_IDWs
  }
   
  if ('ok' %in% intMethods) {
    ## ordinary kriging
    v_mod_ok <- automap::autofitVariogram(elev~1,input_data =
                                            trainingData$sp)
    
    # not happy running with future_lapply!
    interp_OKs <- list()
    for (y in 1:length(paramData.c$ok$run_no)) {
      x <- paramData.c$ok$run_no[y]
      pdata <- paramData.c$ok %>% 
        filter(run_no == x)
      # print(OK_grid[x,])
      st <- Sys.time()
      fit_OK = gstat(formula = elev ~ 1, model = v_mod_ok$var_model, 
                     data = trainingData$sp,
                     nmax=pdata$nmaxVals,
                     nmin=pdata$nminVals)
      interp_OKs[[y]] <- raster::interpolate(testData$ras[[1]], fit_OK)
      tdiff <- Sys.time()-st
      t <- list(val = tdiff,
                unit_chr = units(tdiff))
      intTimes$OK[[y]] <- t
    }
    rasterlist$`Ordinary Kriging` <- interp_OKs
  }
    
  # grass params
  gdb <- paste0(gLoc,Sys.getpid(),'/')
  
  # @@@@ needed to prevent freezing on run 279 ?!? !!!!
  # paramData.c$gspline <- paramData.c$gspline %>%
  #   slice(10:nrow(.))
  # grass-based splines model ----
  # training.sf <- trainingData$sf
  # test.r <- testData$ras
  if ('gspline' %in% intMethods) { 
    # write points for gspline
    st_write(trainingData$sf[,'elev'], 
             paste0('vector/intout_',maskPoly$fid,'_training.gpkg'),
             delete_dsn=T)
    vecLoc <- paste0(getwd(),'/vector/intout_',maskPoly$fid,'_training.gpkg')
    
    system(paste0(
      'grass ',gdb,' --exec g.proj datum=osgb36 -c'
    ))
    
    if (dir.exists(gdb)) unlink(gdb,recursive=T)
    
    system(paste0(
      'grass -c ',gdb,' --exec v.in.ogr input=',vecLoc,' output=points -o --o'
    ))
    
    m <- mask(testData$ras[[1]],trainingData$ras[[1]],inverse=T)
    mLoc <- paste0(getwd(),'/raster/intout_',maskPoly$fid,'_ras_testmask.tif')
    writeRaster(m, mLoc,
                overwrite=T)
    system(paste0(
      'grass ',gdb,' --exec r.in.gdal input=',mLoc,' output=testMask -o --o'
    ))
    system(paste0(
      'grass ',gdb,' --exec g.region raster=testMask'
    ))
    # y <- 1
    interp_GSPLINEs <- list()
    for (y in 1:length(paramData.c$gspline$run_no)) {
      x <- paramData.c$gspline$run_no[y]
      pdata <- paramData.c$gspline %>% 
        filter(run_no == x)
      st <- Sys.time()
      system2('grass',
              paste(shQuote(gdb),
                    shQuote('--exec'),
                    shQuote(paste0(getwd(),'/python/GRASS_vrst.py')),
                    shQuote(pdata$smoothVals),
                    shQuote(pdata$tensionVals),
                    shQuote(pdata$nminVals),
                    shQuote(x),
                    shQuote(pdata$intpol_fid),
                    shQuote(paste0(getwd(),'/raster'))
              ),
              stderr = paste0(getwd(),'/logs/grass_vrst_errout_',
                              pdata$intpol_fid,'.txt')
      )
      print('reading grass output raster...')
      r <- raster(paste0('raster/gspline_int_intfid_',pdata$intpol_fid,
                         '_runnum_',x,'.tif'))
      interp_GSPLINEs[[y]] <- raster::merge(r,trainingData$ras[[1]])
      file.remove(paste0('raster/gspline_int_intfid_',pdata$intpol_fid,
                         '_runnum_',x,'.tif'))
      tdiff <- Sys.time()-st
      t <- list(val = tdiff,
                unit_chr = units(tdiff))
      intTimes$GSPLINE[[y]] <- t
    }
    rasterlist$`GRASS Regularized Splines Tension` <- interp_GSPLINEs
  }
  
  
  if ('gfilter' %in% intMethods) {
    # prepare rasters
    mLoc <- paste0(getwd(),'/raster/intout_',maskPoly$fid,'_ras_training.tif')
    # not sure why rasters merged??
    allRas <- merge(trainingData$ras[[1]],testData$ras[[1]]) # weird?
    # overwritten below:
    allRas <- trainingData$ras[[1]]
    
    writeRaster(allRas, 
                mLoc,
                overwrite=T)
    
    if (dir.exists(gdb)) unlink(gdb,recursive=T)
    
    system(paste0(
      'grass -c ',gdb,' --exec r.in.gdal input=',mLoc,' output=allRas -o --o'
    ))
    system(paste0(
      'grass ',gdb,' --exec g.region raster=allRas'
    ))
    # x <- 1
    # y <- 1
    load(file='cvdev/gfilter.RDS')
    interp_GFILTERs <- list()
    for (y in 1:length(paramData.c$gfilter$run_no)) {
      x <- paramData.c$gfilter$run_no[y]
      
      pdata <- paramData.c$gfilter %>% 
        filter(run_no == x)
      
      # get filter text info
      pdata$filtText <- 
        gfilter.params[which(gfilter.params$filtcomb_number == pdata$filtVals),'col3']
      pdata$filtLength <-
        unlist(lapply(str_split(pdata$filtText,','),length))
      
      pdata$nradVals <- 
          ifelse(pdata$filtLength == 2, 
                 paste0(rep(pdata$radVals,2),collapse=','),
                 as.character(pdata$radVals))
        
      st <- Sys.time()
      system2('grass',
              paste(shQuote(gdb),
                    shQuote('--exec'),
                    shQuote(paste0(getwd(),'/python/GRASS_resampfilter.py')),
                    shQuote(pdata$nradVals),
                    shQuote(pdata$filtText),
                    shQuote(x),
                    shQuote(pdata$intpol_fid),
                    shQuote(paste0(getwd(),'/raster'))
              ),
              stderr = paste0(getwd(),'/logs/grass_resampfilter_errout_',
                              pdata$intpol_fid,'.txt')
      )
      r <- readAll(raster(paste0('raster/gfilter_int_intfid_',pdata$intpol_fid,
                         '_runnum_',x,'.tif')))
      crs(r) <- crs(testData$ras[[1]])
      
      interp_GFILTERs[[y]] <- r
      file.remove(paste0('raster/gfilter_int_intfid_',pdata$intpol_fid,
                         '_runnum_',x,'.tif'))
      tdiff <- Sys.time()-st
      t <- list(val = tdiff,
                unit_chr = units(tdiff))
      intTimes$GFILTER[[y]] <- t
    }
    rasterlist$`GRASS Resampled Filter` <- interp_GFILTERs
  }
  
  # gen list
  rasterlist <- rasterlist %>% 
    map( ~map(.x, .f = function(r) {
      mask(crop(extend(r,testData$ras[[1]]),testData$ras[[1]]),
           testData$ras[[1]])
    }))
  
  # output parameters as melted df
  
  paramsCv <- paramData.c %>% 
    map_df(~reshape2::melt(.x,
                           id.vars=c('run_no','intpol_fid',
                                     'int_method')))
  
  print(paste0('completed pol id...',maskPoly$fid))
  intA <- list(ras = rasterlist,
              params.m = paramsCv,
              intTimes = intTimes)
  
  dat <- compareInt(intRasters=intA,
                    foldedRas=pd$foldA,
                    tiledRas=pd$tiles)
  fout <- paste0(outputDir,'/intdat_',outputTag,'_polfid',pd$pol$fid,'.RDS')
  save(dat,
       file=fout)
  print(paste0('saved file to... ',fout,
         outputDir,'/intdat_',outputTag,'_polfid',pd$pol$fid,'.RDS'))
  return(intTimes)
}

# grass bicubic ----
# tag <- sessionTag
interpolateRasBicubic <- function(pd,cvg,
                                  outputDir = '/media/mal/working_files/quarry/',
                                  testCV=T,
                                  tag) {
  tag <- str_replace(tag,tag,paste0(tag, '_bicubic'))
  trainingData <- pd$foldA$train
  testData <- pd$foldA$all
  maskPoly <-  pd$tiles$pol # if using offset poly
  paramData <- cvg
  
  paramData.c <- 
    paramData %>% 
    map2(.y = names(paramData), 
         .f = function(df, y) {
           if (any(str_detect(names(df),'nmaxVals')) && 
               any(str_detect(names(df),'nminVals'))) {
             df <- df %>% filter(nminVals < nmaxVals) }
           df <- df %>% 
             mutate(run_no = 1:nrow(.),
                    intpol_fid = maskPoly$fid,
                    int_method = y)
           if (testCV) df <- df[1:5,]
           return(df)
         })
  
  
  trainLoc <- paste0(getwd(),'/raster/intout_',maskPoly$fid,'_ras_trainmask.tif')
  writeRaster(trainingData$ras[[1]], trainLoc,
              overwrite=T)
  
  interp_GBICUBICs <- lapply(paramData.c$gbicubic$run_no, function(x) {
    pdata <- paramData.c$gbicubic %>% 
      filter(run_no == x)
    system2('grass',
            paste(
              # shQuote(gdb),
              shQuote('--tmp-location'),
              shQuote('EPSG:27700'),
              shQuote('--exec'),
              shQuote(paste0(getwd(),'/python/GRASS_bspline.py')),
              shQuote(trainLoc),
              shQuote(pdata$stepVals),
              shQuote(pdata$lamVals),
              shQuote(x),
              shQuote(pdata$intpol_fid)
            ),
            stderr = paste0(getwd(),'/logs/grass_bspline_errout_',
                            pdata$intpol_fid,'.txt')
    )
    r <- raster(paste0('raster/gbicubic_int_intfid_',pdata$intpol_fid,
                       '_runnum_',x,'.tif'))
    r.merge <- raster::merge(r,trainingData$ras[[1]])
    file.remove(paste0('raster/gbicubic_int_intfid_',pdata$intpol_fid,
                       '_runnum_',x,'.tif'))
    return(r.merge)
  })
  
  # gen list
  rasterlist <- list(
    "GRASS Bicubic Spline" = interp_GBICUBICs
  ) %>% 
    map( ~map(.x, .f = function(r) {
      mask(crop(extend(r,testData$ras[[1]]),testData$ras[[1]]),
           testData$ras[[1]])
    }))
  
  intA <- list(ras = rasterlist)
  
  
  dat <- compareInt(intRasters=intA,
                    foldedRas=pd$foldA,
                    tiledRas=pd$tiles)
  dat$diff.maps <- NULL
  gc()
  frem <- list.files('raster',pattern=paste0('gbicubic_int_intfid_',pd$tiles$pol$fid),
                     full.names = T)
  file.remove(frem)
  save(dat,
       file=paste0(outputDir,'intdat_',tag,'_polfid',pd$pol$fid,'.RDS'))
  return(rasterlist)
  
  
  
}


# intRasters <- intA
# foldedRas <- foldA
# compareRas <- tiles$b
# maskPoly <- pol

compareInt <- function(intRasters, # list of interpolated rasters
                       foldedRas, # the test/training rasters (raster A, latest)
                       tiledRas # another raster to compare int surfaces with (raster B, earliest)
) { 
  
  # the earlier survey to compare against
  compareRas <- tiledRas$b
  # a polygon to exclude from testErr calculations
  maskPoly <- tiledRas$pol
  
  frtest <- as(foldedRas$test,'Raster')
  
  foldedRas$test <- frtest[[1]]
  
  if (class(compareRas)=='stars') {
    compareRas.r <- as(compareRas, 'Raster')
    compareRas.r <- compareRas.r[[1]]
  }
  
  # this ensures the extent of the comparison raster (older survey)
  # is the same as the interpolated surfaces. Errors can occur when
  # earlier survey did not quite cover the same cells as the later one.
  compareRas.r <- mask(crop(extend(compareRas.r,foldedRas$all$ras[[1]]),
                 foldedRas$all$ras[[1]]),
            foldedRas$all$ras[[1]])
  
  # training error is difference between original training data and modeled training data
  # test error is difference between original test data and modeled test data
  
  # test data below is essentially the original raster 
  # ep <- pd$pol
  # fr <- foldedRas
  # cr <- compareRas.r
  compareEach <- function(raslist,fr,cr,ep=NULL) {
    
    # raslist <- intRasters$ras$`GRASS Regularized Splines Tension`
    # interpolated <- raslist[[1]]
    compFunction <- function(interpolated) {
      
      trainingErr.r <- mask(interpolated,fr$train$ras[[1]]) - 
        fr$train$ras[[1]]
      
      testErr.r <- mask(interpolated,fr$test) - fr$test
      compareDiff <- interpolated - cr
      
      out <- list(
        trainingErr.r = trainingErr.r,
        testErr.r = testErr.r,
        compareDiff = compareDiff)
      
      if (!is.null(ep)) {
        ep.r <- fasterize::fasterize(ep, fr$test)
        out$testErr.ex.r <- mask(interpolated,fr$test) - 
          mask(fr$test,ep.r,inverse=T)
        out$testErr.inc.r <- mask(interpolated,fr$test) - 
          mask(fr$test,ep.r)
        out$compareDiff.inc.r <- interpolated - mask(cr,ep.r)
        out$compareDiff.ex.r <- interpolated - mask(cr,ep.r,inverse=T)
      }
      
      return(out)
    }
    
    updatedRas <- raslist %>% map(~compFunction(.x))
    # raslist$ras.comp <- updatedRas
    # return(raslist)
  }
  
  diffMaps <- intRasters$ras %>% 
    map(~compareEach(.x, 
                     fr = foldedRas,
                     cr = compareRas.r,
                     ep = maskPoly))
  
  # diffMaps$`Nearest Neighbor`$compare$compareDiff
  # plot(stack(diffMaps$`Triangular Irregular Surface`))
  # r <- diffMaps$`Ordinary Kriging`$trainingErr.r
  calcRMSEfromRas <- function(r) sqrt(cellStats(r^2,mean))
  
  diffRMSEs <- lapply(names(diffMaps), function(x) {
    diffRMSEs <- diffMaps[[x]] %>% map_df(map, calcRMSEfromRas) %>% 
      mutate(int_method = x,
             run_no = 1:nrow(.))
  }) %>% bind_rows()
  
  # add pol fid if pol provided
  if (!is.null(maskPoly)) diffRMSEs <- diffRMSEs %>% 
    mutate(intpol_fid = maskPoly$fid)
  
  return(list(orig.maps = intRasters,
              # diff.maps = diffMaps,
              tiles = tiledRas,
              rmses = diffRMSEs))
}

visData <- function(dat) {
  rmses.df <- dat$rmses
  
  # first plots of CV data (if it exists)
  if (nrow(dat$orig.maps$`Nearest Neighbor`$param) > 1) {
    params <- dat$orig.maps %>% map(pluck,'param')
    params[unlist(lapply(params,function(x) length(x) == 0))] <- NULL
    flatten(params)
    params %>% map_df(flatten)
    rmses.df
  }
  
  dat$orig.maps %>% map_df('param') 
  
  rmsesByMethod <- split(dat$rmses, dat$rmses$int_method)
  rmsesByMethod %>% map(left_join, .x, )
  
  
}

# cross validation check
# ras <- interpRas[[3]]$rasterListClip$`Nearest Neighbor`$surf$orig
# tensionVals <- 10
# smoothVals <-  0.5
# npminVals <-  50
# ras <- rasPrep$hole$hole.st.elev
crossValidateSplines <- function(ras, tensionVals = seq(0.01,0.1,by=0.01),
                                 smoothVals = seq(10,30,by=5),
                                 npminVals = seq(60,140,by=20),
                                 tag='sys.time') {
  
  if (tag == 'sys.time') {
    tag <- paste0('t_',format(Sys.time(),'%d_%m_%Y_%H_%M_%S')) }
  
  rasSf <- ras %>% 
    st_as_sf(points=TRUE,merge=FALSE) %>% 
    st_centroid() %>% 
    mutate(cat = 1:nrow(.)) %>% 
    dplyr::rename(elev = 1)
  
  # rasSf <- rasSf %>% # test
  #   head(1000)
  
  st_write(rasSf,'vector/op_points.gpkg',delete_dsn=T)
  # ras <- rasterFromXYZ(cbind(rasSf %>% st_coordinates(),rasSf$layer)) # test
  
  maskLoc <- 'raster/op_raster.tif'
  # writeRaster(ras,maskLoc,overwrite=T) # test
  writeRaster(as(ras,"Raster"),
              file=maskLoc,
              overwrite=T)
  
  pointsLoc <- paste0(getwd(),'/vector/op_points.gpkg')
  
  # optimise these parameters
  
  # starter params
  # tensionVals <- seq(10,200,by=10)
  # smoothVals <- seq(0.1,1,by=0.1)
  # npminVals <- seq(200,400,by=100)
  # best run with the above was
  # cvdevGrid[116,]
  #         er_mean     er_mean_abs run smoothVal tensionVal npminVal
  # 116 -1.7246e-05   0.0159536   204       0.4         10      300
  
  pGrid <- expand.grid(smoothVal = smoothVals,
                       tensionVal = tensionVals,
                       npminVal = npminVals) %>% 
    mutate(run = 1:nrow(.))
  
  plan(multisession, gc = TRUE, workers = 4) ## Run in parallel on local computer
  # x <- 1
  
  future.apply::future_lapply(1:nrow(pGrid), function(x) {
    
    system2('grass',
            paste(shQuote('--tmp-location'),
                  shQuote('EPSG:27700'),
                  shQuote('--exec'),
                  shQuote('/home/barneyharris/projects/quarry/python/GRASSoptimise.py'),
                  shQuote(pointsLoc),
                  shQuote(pGrid[x,'smoothVal'] ),
                  shQuote(pGrid[x,'tensionVal']),
                  shQuote(pGrid[x,'npminVal']),
                  shQuote(x),
                  shQuote(maskLoc),
                  shQuote(tag)
            ),
            stderr = paste0(getwd(),'/logs/grassOptimise_errout.txt')
    )
  })
  
  # x <- list.files('cvdev',full.names = T)[1]
  # x <- "cvdev/cvdev_1_t_16_12_2020_13_30_44.txt"
  cvdevs <- do.call(rbind,lapply(list.files('cvdev',
                                            pattern=tag,
                                            full.names = T), function(x) {
                                              runNum <- as.numeric(gsub(".*?([0-9]+).*", "\\1", x))
                                              d <- read.csv(x, stringsAsFactors = F, header=F)
                                              vals <- lapply(d$V1, str_split_fixed, '=', 2) %>% map(2)
                                              d$vals <- as.numeric(unlist(vals))
                                              # data.frame(var = c('mean','mean_abs'),
                                              #            val = c(d[8,2],d[9,2]),
                                              #            run = runNum)
                                              data.frame(mae = d[9,2],
                                                         rmse = sqrt(d[24,2]),
                                                         run = runNum)
                                            }))
  
  cvdevGrid <- left_join(cvdevs,pGrid)
  
  
  return(cvdevGrid)
}

# functions for 3_analysing results -----

# function for plotting the location data, showing earlier survey,
# later survey and difference between surveys
# pdat <- prepData.c[[1]]
# TO DO
# add location X / Y to plot title
# add scale bar to plots
# reverse subtraction of A and B surveys
plotPreppedData <- function(pdat,silent=F) {
  
  # labels
  my_params <- list(scale_fill_viridis(na.value="transparent"),
    coord_sf(datum = sf::st_crs(27700)))
  my_theme <- theme(axis.text=element_blank(),
                    axis.title=element_blank(),
                    axis.ticks = element_blank())
  my_g <- theme_grey(base_size = 8) 
                         
  
  if (!silent) print(pdat$pol)
  df <- pdat$tiles$b %>% as.data.frame() %>% 
    dplyr::rename(value = 3)
  
  bMap <- ggplot() +  
    geom_raster(data=df, 
                aes(x=x, y=y, fill=value)) + 
  labs(title = paste0('Site: ',sprintf("%03d", unique(pdat$pol$fid))),
       subtitle = paste0('Survey date: ',unique(pdat$pol$b))) + 
    my_params + my_g + my_theme
    
  
  df <- pdat$tiles$a %>% as.data.frame() %>% 
    dplyr::rename(value = 3)
  
  aMap <- ggplot() +  
    geom_raster(data=df, 
                aes(x=x, y=y, fill=value)) + 
    labs(title = paste0('Site: ',sprintf("%03d", unique(pdat$pol$fid))),
         subtitle = paste0('Survey date: ',unique(pdat$pol$a))) + 
    my_params + my_g + my_theme
  
  if (!is.null(pdat$abDiff$abdiff_st)) {
    df <- pdat$abDiff$abdiff_st %>% as.data.frame() %>% 
      dplyr::rename(value = 3)
    
    diffMap <- ggplot() +  
      geom_raster(data=df, 
                  aes(x=x, y=y, fill=value)) + 
      labs(title = paste0('Site: ',sprintf("%03d", unique(pdat$pol$fid))),
           subtitle = paste0('RMSE: ',
                             as.character(round(as.numeric(unique(pdat$abDiff$abdiff_rmse)),4))
                             )) + 
      my_params + my_g + my_theme
    
    multi.page <- ggpubr::ggarrange(plotlist = list(bMap,aMap,diffMap), 
                                    nrow = 1, ncol = 3)
    return(multi.page)
  } else return(NULL)
  
}


# filePattern="intdat_prepData_bestOffSets_buffpol10_smper0_all"
# filePattern="intdat_prepData_bestOffSets_buffpol10_smper50_all"
# filePattern <- 'prepData_unmod_locs_NYSE_2013_2009_1to50_smpper0_all'
# plot2pdf <- F
# descPlots <- F
# bestPlots = F
# oldSave = T
# dirLoc = '/media/mal/working_files/quarry'
analyseDat <- function(dirLoc = '/media/mal/working_files/quarry',
                       filePattern,
                       descPlots = T, # error plots, var plots etc
                       bestPlots = T, # calc best run & plot
                       plot2pdf = T,
                       oldSave = F) { 
  # loading function
  loadRData <- function(fileName) {
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  
  # loads pre-saved RDS files from dirLoc locations ...
  fs <- list.files(dirLoc,full.names = T,
                   pattern=paste0(filePattern,'_pol*'))
  # fs <- list.files(paste0(dirLoc,'/backup'),full.names = T,
  #                  pattern=paste0(filePattern,'_pol*'))
  # 
  
  print('loading results...')
  if (oldSave) {
    rasData <- lapply(fs, loadRData) } else {
      rasData <- lapply(fs, readRDS) }
  print('done')
  
  
  # rdFilter <- rasData
  # replace old resamp data rasters
  # for (p in 1:length(rasData)) {
  #   rasData[[p]]$orig.maps$ras$`GRASS Resampled Filter` <-
  #     rdFilter[[p]]$orig.maps$ras$`GRASS Resampled Filter`
  # }
  # # replace old resamp rmses tables
  # for (p in 1:length(rasData)) {
  #   newrmse <- rasData[[p]]$rmses %>%
  #     dplyr::filter(int_method != "GRASS Resampled Filter")
  #   newrmse <- rbind(newrmse,rdFilter[[p]]$rmses)
  #   rasData[[p]]$rmses <- newrmse
  # }
  # 
  # for (x in 1:length(fs)) {
  #   print(x)
  #   s <- rasData[[x]]
  #   save(s,file=fs[x])
  # }
  
  # name correlations
  nCor <- data.frame(long_name = c("Nearest Neighbor",
                                   "Inverse Distance Weighted",
                                   "Ordinary Kriging",
                                   "GRASS Regularized Splines Tension",
                                   "GRASS Bicubic Spline",
                                   "GRASS Resampled Filter"),
                     short_name = c('nn','idw','ok','gspline',
                                    'gbicubic','gfilter'))
  
  # extract rmses for each pol
  rmsesAll <- rasData %>% 
    map_df( ~.x$rmses)
  
  # can check specific runs...
  # check <- rmsesAll %>%
  #   dplyr::filter(int_method == "GRASS Resampled Filter" &
  #                   intpol_fid == 6 &
  #                   run_no == 231)
  
  # processes resampled maps to remove runs with holes etc.
  
  if (any(rmsesAll$int_method %in% "GRASS Resampled Filter")) {
    # # for resamp filter, remove any ras with holes in
    rd <- rasData %>% map(.f = function(p) {
      origNas <- as.numeric(table(is.na(p$tiles$a[[1]][]))[2])
      p$orig.maps$ras$`GRASS Resampled Filter` <-
        p$orig.maps$ras$`GRASS Resampled Filter` %>%
        map(.f = function(r) {
          intNas <- as.numeric(table(is.na(r[]))[2])
          if (intNas == origNas) r
        })
      # %>% purrr::compact()
      p
    })
    
    # get df of name and run id of valid maps
    validRuns <- rd %>% map(.f = function(p) {
      p$orig.maps$ras$`GRASS Resampled Filter` %>%
        purrr::compact() %>%
        map_chr(.f = function(x) names(x))
    }) %>% unlist() %>% data.frame(rasname=.) %>%
      tidyr::separate(rasname, as.character(1:6), "_") %>%
      dplyr::select(4,6) %>%
      dplyr::rename('intpol_fid' = '4',
                    'run_no' = '6') %>%
      dplyr::mutate_if(is.character,as.numeric) %>%
      dplyr::mutate(int_method = 'GRASS Resampled Filter')
  
    validResamp <- rmsesAll %>% 
      dplyr::filter(int_method == 'GRASS Resampled Filter') %>%
      dplyr::semi_join(validRuns)
    
    rmsesAll <- rmsesAll %>% 
      dplyr::filter(int_method != 'GRASS Resampled Filter') %>% 
      rbind(validResamp)
  }

  # group vars
  gVars = c('run_no','intpol_fid','int_method')
  # err vars - variables to be plotted
  errVars <- setdiff(names(rmsesAll),gVars)
  
  # if all testErr.ex.r is NaN (= no noise added) then replace with 0
  rmsesAll[is.nan(rmsesAll$testErr.ex.r),'testErr.ex.r'] <- 
    0
                        
  # filter out rows relating to failed interpolations
  rmsesAll <- rmsesAll %>% 
    dplyr::filter(across(all_of(c(setdiff(errVars,'testErr.ex.r'))),
                         ~ !is.nan(.)))
  
  # error box plots ---- 
  # these plots show the range of error values associated by all interpolations
  # grouped by int method and site ID
  if (descPlots) {
  print('plotting error plots...')
    errPlots <- errVars %>% map(.f = function(v) {
      rmsesAll2 <- rmsesAll
      rmsesAll2$intpol_fid <- 
        as.character(rmsesAll2$intpol_fid)
      
      g <- ggplot(rmsesAll2) + 
        geom_boxplot(aes_string(x = 'intpol_fid',
                                y = v)) + 
        facet_wrap(~int_method,scales='free_y') + 
        theme(axis.text.x = element_text(angle = 45, hjust=1))
      
    }) %>% setNames(errVars)
    
    if (plot2pdf) {
      print('exporting pdf...')
      multi.page <- ggpubr::ggarrange(plotlist = errPlots, 
                                      nrow = 1, ncol = 1) 
      ggpubr::ggexport(multi.page, 
                       filename =
                         paste0('plots/',filePattern,'_errorPlots.pdf'),
                       height = 8.25, width = 11.75)
    }
  } # /if (descPlots)
  
  errCoVars <- 
    data.frame(c())
  if (descPlots) {
    errCovarPlot <- errVars %>% map(.f = function(v) {
      errVars
    })
  }
  
  # variable plots ----
  allData <- rasData %>% 
    map_df( ~.x$orig.maps$params.m) %>% 
    left_join(nCor, by=c('int_method' = 'short_name')) %>% 
    mutate(int_method = long_name) %>% 
    dplyr::select(-long_name) %>% 
    left_join(rmsesAll, by=c('run_no','intpol_fid','int_method'))
  
  # plots showing how cross validation parameters result in differing error
  # values. Parameters with only a single value across all interpolations
  # e.g. nMinVals will show as blank plots
  # erv <- errVars[2]
  # i <- 'GRASS Bicubic Spline'
  if (descPlots) {
  print('plotting var plots...')
    varPlots <- errVars %>% map(.f = function(erv) {
      
      varPlots <- unique(allData$int_method) %>% 
        map(.f = function(i) {
          
          df <- allData %>% dplyr::filter(int_method == i)
          
          g <- ggplot(df) + 
            stat_smooth(aes(x = value,
                            y = get(erv)),
                        method = 'loess') +
            facet_wrap_paginate(intpol_fid ~ variable,
                                scales = 'free',
                                ncol = length(unique(df$variable)), # n variables
                                nrow = 4, # n sites
                                page = 1) + 
            ggtitle(i) + 
            ylab(erv) + 
            xlab('Parameter value') + 
            theme_grey(base_size=8)
          return(g)
        }) %>% setNames(unique(allData$int_method))
    }) %>% setNames(errVars) # / varPlots
  
    exportVars <- c('testErr.ex.r','testErr.inc.r',
                    'compareDiff.ex.r',
                    'compareDiff.inc.r')
    
    if (plot2pdf) {
      exportVars %>% map(.f = function(ex) {
        multi.page <- ggpubr::ggarrange(plotlist = varPlots[[ex]], 
                                        nrow = 1,ncol=1) 
        
        ggpubr::ggexport(multi.page, 
                         filename =
                           paste0('plots/',filePattern,'_',ex,'_varPlots.pdf'),
                         height = 8.25, width = 11.75)
        })
    }
  } #  / if (descPlots)
  
  # plots showing relaitionship between test error values 
  # and earlier surface values
  
  # !! plots may well be empty if no random noise introduced in preparation of
  # data using prepareData(smpper = N) argument
  
  # cross plots for testErr.ex.r ----
  rmsesAll.m <- rmsesAll %>% 
    reshape2::melt(id.vars = c(gVars, 'testErr.ex.r'))
  
  # vn <- 1
  if (descPlots) {
  print('plotting cross plots ex.r...')
    crossPlots.testErr.ex.r <-  1:length(unique(rmsesAll.m$variable)) %>% 
      map(.f = function(vn) {
        df <- rmsesAll.m %>% 
          dplyr::filter(variable == unique(rmsesAll.m$variable)[vn]) %>% 
          mutate_if(is.factor,as.character)
        
        # i <- allData$int_method[1]
        p <- unique(allData$int_method) %>% 
          map(.f = function(i) {
            ggplot(df %>% dplyr::filter(int_method == i)) +
              stat_smooth(aes(x = testErr.ex.r,
                              y =  value),
                          method = 'loess') + 
              facet_wrap(~ intpol_fid,
                         scales = 'free',
                         # ncol = 2, # n variables
                         nrow = 3 # n sites
              ) + 
              ggtitle(paste(i,unique(rmsesAll.m$variable)[vn],
                            collapse=': ')) +
              ylab(unique(rmsesAll.m$variable)[vn])
          }) %>% setNames(unique(allData$int_method))
      }) %>% setNames(unique(rmsesAll.m$variable)) # / crossPlots
    
    if (plot2pdf) {
      # export pdf
      names(crossPlots.testErr.ex.r) %>% 
        map(.f = function(v) {
          multi.page <- ggpubr::ggarrange(plotlist = crossPlots[[v]], 
                                          nrow = 1, ncol = 1) 
          ggpubr::ggexport(multi.page, 
                           filename =
                             paste0('plots/',filePattern,'_',v,
                                    '_crossPlots.testErr.ex.r.pdf'),
                           height = 8.25, width = 11.75)
        })
    }
  } # / if (descPlots)
  
  # cross plots for testErr.inc.r ----
  rmsesAll.m <- rmsesAll %>% 
    reshape2::melt(id.vars = c(gVars, 'testErr.inc.r'))
  
  # vn <- 1
  if (descPlots) {
    print('plotting error plots inc.r...')
    crossPlots.testErr.inc.r <-  1:length(unique(rmsesAll.m$variable)) %>% 
      map(.f = function(vn) {
        df <- rmsesAll.m %>% 
          dplyr::filter(variable == unique(rmsesAll.m$variable)[vn]) %>% 
          mutate_if(is.factor,as.character)
        
        # i <- allData$int_method[1]
        p <- unique(allData$int_method) %>% 
          map(.f = function(i) {
            ggplot(df %>% dplyr::filter(int_method == i)) +
              stat_smooth(aes(x = testErr.inc.r,
                              y =  value),
                          method = 'loess') + 
              facet_wrap(~ intpol_fid,
                         scales = 'free',
                         # ncol = 2, # n variables
                         nrow = 3 # n sites
              ) + 
              ggtitle(paste(i,unique(rmsesAll.m$variable)[vn],
                            collapse=': ')) +
              ylab(unique(rmsesAll.m$variable)[vn])
          }) %>% setNames(unique(allData$int_method))
      }) %>% setNames(unique(rmsesAll.m$variable)) # / crossPlots
    
    if (plot2pdf) {
      # export pdf
      names(crossPlots.testErr.inc.r) %>% 
        map(.f = function(v) {
          multi.page <- ggpubr::ggarrange(plotlist = crossPlots[[v]], 
                                          nrow = 1, ncol = 1) 
          ggpubr::ggexport(multi.page, 
                           filename =
                             paste0('plots/',filePattern,'_',v,
                                    '_crossPlots.testErr.inc.r.pdf'),
                           height = 8.25, width = 11.75)
        })
    }
  } # / if (descPlots)
  
  # best runs routine ----
  # x <- errVars[5]
  # nRuns <- 2
  print('calculating best runs...')
  getNbestRuns <- function(errVars, # which error values to process
                           nRuns=1  # how many best runs to return
                           ) {
      errVars %>% map_df(.f = function(x) {
        # calc quants
        # find best runs per error value - this function groups data by site &
        # int method and splits each error metric into quantiles. Returns the
        # value below which represents the best 10% of runs
        print(paste0('finding runs with lowest errors for...',x))
        rmsesQ <- rmsesAll %>% 
          dplyr::group_by(intpol_fid,int_method) %>% 
          dplyr::summarise(quantile = scales::percent(c(0.1)),
                    qv = 
                      quantile(get(x), c(0.1)))
        # run ID and exact error value of run(s) are returned, number of 
        # runs can be specified with 'nRuns'
        a <- rmsesAll %>% left_join(rmsesQ) %>% 
          dplyr::group_by(intpol_fid,int_method) %>% 
          dplyr::filter(get(x) <= qv) %>% 
          { if (!is.null(nRuns)) slice_min(., order_by=get(x), n=all_of(nRuns)) else .} %>% 
          dplyr::select(-setdiff(errVars,x)) %>% 
          dplyr::rename('error_value' = all_of(x)) %>% 
          mutate(error_value = as.numeric(error_value),
                 error_var = all_of(x) )
        }) %>%
      # dump traininErr.r as most models have zero error for this
      dplyr::filter(error_var != 'trainingErr.r') %>% 
      left_join(allData[,c(gVars,'variable','value')] %>% 
                  filter(variable == 'nmaxVals')) %>% 
      { if (!is.null(nRuns))  group_by(.,intpol_fid,int_method, error_var) %>% 
          slice_min(order_by = 'nmaxVals') else . } %>% 
      # select run with lowest nmax (reduces computation time)
      mutate(bid = 1:nrow(.)) %>% 
      dplyr::select(-c('variable', 'value'))
  }
  # find best single run per pol, int method and erro var
  bestRuns <- getNbestRuns(errVars, nRuns=1)
  
  # isoloate parameters for best runs, returned in long format
  bestRuns.params <- bestRuns %>% 
    left_join(allData,
              by = c("int_method", "run_no", "intpol_fid")) %>% 
    dplyr::select(c(all_of(gVars),'bid','variable','value'))
  
  # find top quartile of runs and calculate range of parameter values
  errVars.priority <- c("testErr.ex.r",
                        "testErr.inc.r",
                        "compareDiff.inc.r")
  # same as above but for limited number of error metrics
  topQRuns <- getNbestRuns(errVars.priority, 
                           nRuns=1)
  # same as above...
  topQRuns.params <- topQRuns %>% 
    left_join(allData,
              by = c("int_method", "run_no", "intpol_fid")) %>% 
    dplyr::select(c(all_of(gVars),'bid','variable','value'))
  
  # construct a parameter range using min and max for 
  # parameters associated with best runs for selected
  # error metrics
  topQRuns.params.range <- topQRuns.params %>% 
    na.omit() %>% 
    group_by(int_method, variable) %>% 
    dplyr::summarise(min_value = min(value),
                     max_value = max(value))
  # split into a list based on the int method and make data wider to
  # compare parameters across different runs
  topQRuns.params.intlist <- topQRuns.params %>% 
    split(f=.$int_method) %>% 
    map(~tidyr::pivot_wider(.x, 
                    names_from = variable,
                    values_from = value)) %>% 
    map(~as.data.frame(left_join(.x, rmsesAll, by=c('intpol_fid','run_no',
                                      'int_method'))))
  
  # error bar plots -----
  errBars <- ggplot(bestRuns) + 
    geom_col(aes(x = as.character(intpol_fid),
                 y = error_value,
                 fill = int_method),
             position='dodge'
    ) + 
    facet_wrap(~ error_var) + 
    theme(legend.position = 'bottom')
  
  # error boxplots plots -----
  bestRuns.f <- bestRuns %>% 
    dplyr::filter(error_var %in% 
                    c('compareDiff.inc.r',
                      'testErr.inc.r'))
  errBox <- ggplot(bestRuns.f) + 
    geom_boxplot(aes(x = int_method,
                 y = error_value)) + 
    facet_wrap(~ error_var, scales='free_y') + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = -45, hjust = 0))
  
  errBox
  # extract best rasters and calculate difference maps routine -----
  
  # begin by converting raster data names to reflect real polygon fids
  names(rasData) <- rasData %>% 
    map(~unique(.x$rmses$intpol_fid)) %>% 
    unlist(.)
  
  # extract rasters to flat list and generate difference maps using rasters
  # identified through best run IDs
  # x <- 1
  # as(rasData$`2`$tiles$a,'Raster')
  # plot(rasData$`2`$orig.data$ras$layer.1)
  bestRas <- 1:nrow(bestRuns) %>% 
    map(.f = function(x) {
      df <- bestRuns[x,]
      intras <- rasData[[as.character(df$intpol_fid)]]$orig.maps$ras[[df$int_method]][[df$run_no]]
      diffras_a <- as(rasData[[as.character(df$intpol_fid)]]$tiles$a,'Raster') - intras
      diffras_b <- as(rasData[[as.character(df$intpol_fid)]]$tiles$b,'Raster') - intras
      list(intras = intras,
           diff_a = diffras_a,
           diff_b = diffras_b)
    }) %>% setNames(as.character(1:nrow(bestRuns)))
  
  # split list into logical groupings & # extract original interpolated rasters, 
  # by map with lowest error per error measure
  bestRas.l <- split.data.frame(bestRuns, bestRuns$error_var) %>% 
    map(~split(.x, .x$intpol_fid)) %>% map(~map(.x, .f = function(b) {
      list(ras = bestRas[b$bid],
           params = b %>% 
             left_join(bestRuns.params %>% dplyr::select(-run_no) ))
    }) 
    # %>% setNames(b$int_method)
    )
  
  # prepare polygons for plotting
  allPol <- rasData %>% map(~.x$tiles$pol) %>% 
    bind_rows()
  
  # plot best ras
  # bestRas.l <- results$prepData_bestOffSets_buffpol10_smper0_all$bestRas.l
  # errRas <- bestRas.l[[1]]
  # fidRas <- errRas$`4`
  
  if (bestPlots) {
    print('plotting best raster plots...')
    bestRasPlots <- bestRas.l %>% map(.f= function(errRas) {
      map(errRas, .f = function(fidRas) {
        
        allRasdf <- map2_df(fidRas$ras, .y = names(fidRas$ras), .f = function(r, rn) {
          df <- as.data.frame(rasterToPoints(r$intras)) %>% 
            dplyr::rename(elev = 3) %>% 
            mutate(bid = as.numeric(rn),
                   type = 'int') %>% 
            left_join(bestRuns, by='bid')
          
          df2 <- as.data.frame(rasterToPoints(r$diff_a)) %>% 
            dplyr::rename(elev = 3) %>% 
            mutate(bid = as.numeric(rn),
                   type = 'diff_a') %>% 
            left_join(bestRuns, by='bid')
          
          df3 <- as.data.frame(rasterToPoints(r$diff_b)) %>% 
            dplyr::rename(elev = 3) %>% 
            mutate(bid = as.numeric(rn),
                   type = 'diff_b') %>% 
            left_join(bestRuns, by='bid')
          
          rbind(df,df2,df3)
        })
        # parse parameters and RMSEs
        tmp <- fidRas$params %>%
          dplyr::group_by(int_method, run_no, bid, error_value) %>% 
          dplyr::filter(n() > 1) %>% 
          dplyr::summarise(param_line = paste0(variable,': ',value,collapse='\n')) %>% 
          mutate(param_all = paste0(param_line, '\nRMSE: ',round(error_value,4))) %>% 
          dplyr::select(-param_line)
        
        df.l <- fidRas$params %>%
          dplyr::group_by(int_method, run_no, bid, error_value) %>% 
          dplyr::filter(n() == 1) %>% 
          mutate(param_all = paste0('RMSE: ',round(error_value,4))) %>% 
          dplyr::select(all_of(names(tmp))) %>% 
          bind_rows(tmp)
        
        plotPol <- allPol %>% dplyr::filter(fid == unique(fidRas$params$intpol_fid))
        
        elevPlots <- ggplot() +  
          geom_raster(data=allRasdf %>% filter(type == 'int'), 
                      aes(x=x, y=y, fill=elev)) + 
          geom_sf(data = plotPol, fill=NA,color='black',
                  linetype = "dashed", size = 0.5) +
          scale_fill_viridis(oob=scales::squish) +
          coord_sf(datum = sf::st_crs(27700)) +
          theme(legend.position="none") + 
          facet_wrap(~ int_method) + 
          labs(title = paste0('Interpolated surfaces for site: ',unique(allRasdf$intpol_fid)),
               subtitle = paste0('Lowest error maps for ',unique(allRasdf$error_var))
          )
        
        diffPlots <- c('diff_a','diff_b') %>% 
          map(.f = function(difftype) {
            ggplot() +  
              geom_raster(data=allRasdf %>% filter(type == all_of(difftype)), 
                          aes(x=x, y=y, fill=elev)) + 
              geom_sf(data = plotPol, fill=NA,color='black',
                      linetype = "dashed", size = 0.5) + 
              scale_fill_viridis(oob=scales::squish) +
              coord_sf(datum = sf::st_crs(27700)) + 
              theme(legend.position="bottom",
                    legend.key.width = unit(10,'mm')) + 
              facet_wrap(~ int_method) + 
              labs(title = paste0('Difference maps for site: ',unique(allRasdf$intpol_fid)),
                   subtitle = paste0('Lowest error maps for ',unique(allRasdf$error_var)),
                   fill = 'Difference (m)'
              )
          }) %>% setNames(c('diff_a','diff_b'))
        
        gl <- list(elev = elevPlots,
             diff_a = diffPlots$diff_a,
             diff_b = diffPlots$diff_b)
        
        # check width of raster to make sure plot is wide enough for labels
        minWidth <- 200
        rasWidth <- (max(allRasdf$x) - min(allRasdf$x))
        if (rasWidth < minWidth) {
          
          gl <- lapply(gl, function(g) {
            newXmin <- min(allRasdf$x)-((minWidth-rasWidth))
            newXmax <- max(allRasdf$x)+((minWidth-rasWidth))
            g + lims(x = c(newXmin,newXmax)) + 
              geom_label(data=df.l, label.size=0.1, size=2.75,
                         aes(x=newXmin,
                             y=min(allRasdf$y),
                             label=param_all),
                         nudge_x=30,
                         nudge_y=30) 
          })
        } else {
          gl <- lapply(gl, function(g) {
            g + geom_label(data=df.l, label.size=0.1, size=2.75,
                           aes(x=min(allRasdf$x),
                               y=min(allRasdf$y),
                               label=param_all),
                           nudge_x=(max(allRasdf$x) - min(allRasdf$x))*0.2,
                           nudge_y=(max(allRasdf$y) - min(allRasdf$y))*0.2
                     )
            }) 
        } # else
        return(gl)
      })
    })
    
    # export to pdf
    if (plot2pdf) {
      names(bestRasPlots) %>% 
        map(.f = function(erv) {
          multi.page <- ggpubr::ggarrange(plotlist = flatten(bestRasPlots[[erv]]), 
                                          nrow = 1, ncol = 1) 
          ggpubr::ggexport(multi.page, 
                           filename =
                             paste0('plots/',filePattern,'_',erv,
                                    '_bestRasPlots.pdf'),
                           height = 8.25, width = 11.75)
        })
    }
  } # / if (bestPlots)
  
  plots <- list()
  
  if (descPlots) {
    plots$errPlots <- errPlots
    plots$varPlots <- varPlots
    plots$crossPlots.testErr.inc.r <- crossPlots.testErr.inc.r
    plots$crossPlots.testErr.ex.r <- crossPlots.testErr.ex.r
  }
  
  if (bestPlots) {
    plots$bestRasPlots <- bestRasPlots
  }
  
  # saveRDS(plots$bestRasPlots,
  #         file=paste0(filePattern,'bestRas_plots.RDS'))
  
  l <- list(bestRas.l = bestRas.l,
            bestRuns = bestRuns,
            bestRuns.params = bestRuns.params,
            topQRuns.params.intlist = topQRuns.params.intlist,
            allData = allData,
            filePattern = filePattern)
  
  # if plots haven't been sent to PDF then add them to R object
  # (this creates a large R object size)
  if (!plot2pdf) l$plots <- plots
  
  return(l)
}

# takes input from analyseDat()
# sessionResults <- results$prepData_bestOffSets_buffpol10_smper50_all
# sessionResults <- results$prepData_bestOffSets_buffpol10_smper0_all
# results <-
#   readRDS('data/results_prepData_unmod_locs_NYSE_2013_2009_1to50_smpper0_all.RDS')
# sessionResults <- results[[1]]
exploreFurther <- function(sessionResults) {
  # examine the best interpolation method per site
  
  # use file pattern to load prepared data to get ab diffs
  # res <- paste0('^',str_match(sessionResults$filePattern, 
  #                             "intdat_\\s*(.*?)\\s*_all")[,2],
  #               '.RDS')
  res <- paste0('^',str_replace(sessionResults$filePattern,
                                '_all',''),'.RDS')
  prepData <- readRDS(list.files('data',pattern=res,
                                 full.names = T))
  # add names to prepData if missing
  if (is.null(names(prepData))) {
    n <- prepData %>% map(~unique(.x$pol$fid)) %>% unlist()
    names(prepData) <- n
  }
  
  # analyse A surveys in terms of sloe, tpi etc.
  # n <- 1
  terrainA <- names(prepData) %>% map_df(.f = function(n) {
    r.slope <- terrain(prepData[[n]]$foldA$train$ras,'slope')
    r.tri <- terrain(prepData[[n]]$foldA$train$ras,'TRI')
    r.rough <- terrain(prepData[[n]]$foldA$train$ras,'roughness')
    statVar <- c('mean','max','sd')
    # statv <- 'range'
    rstats <- statVar %>% 
      map(.f = function(statv) {
      df <- c(r.slope,r.tri,r.rough) %>% 
        map_dbl(~cellStats(.x,statv)) %>%
        data.frame() %>% 
        setNames('var_value') %>% 
        mutate(var_name = paste0(statv,'_',c('slope','tri','roughness')),
               `Site ID` = as.numeric(n))
        
      })
  }) %>% 
    tidyr::pivot_wider(names_from = c(var_name),
                       values_from = c(var_value)) 
  
  abDiffs <- names(prepData) %>% map_df(.f = function(n) {
    data.frame(`Site ID` = as.numeric(n),
               `Earlier survey and later survey (polygon): RMSE`= 
                 prepData[[n]]$abDiff$abdiff_rmse.inc.r,
               check.names=F)
  })
  
  # investigate best runs
  
  
  bestInts.list <- sessionResults$bestRuns %>% 
    dplyr::group_by(intpol_fid, error_var) %>% 
    slice_min(order_by = error_value) %>% 
    group_by(error_var) %>% 
    group_split()
  names(bestInts.list) <- unique(sessionResults$bestRuns$error_var)
  
  extractRas <- function(b) {
    errVar <- unique(bestInt$error_var)
    # polras <- results[[1]]$bestRas.l[[errVar]][[1]]
    ras <- names(results[[1]]$bestRas.l[[errVar]]) %>% 
      map(.f = function(polfid) {
        bid <- b %>% dplyr::filter(intpol_fid == polfid) %>% 
          dplyr::select(bid) %>% as.character()
        results[[1]]$bestRas.l[[errVar]][[polfid]]$ras[[bid]]
      }) %>% setNames(names(results[[1]]$bestRas.l[[errVar]]))
  }
  
  bestInt.l.r <- bestInts.list %>% map(extractRas)
  
  # process best run data into table
  
  # first collapse int parameters onto single line
  bestRuns.params.sum <- sessionResults$bestRuns.params %>% 
    dplyr::mutate(params_conc = paste0(variable,' = ',value)) %>% 
    dplyr::group_by(bid) %>% 
    dplyr::summarise(params_conc = 
                       paste0(int_method, ' with ',
                              paste0(params_conc,collapse='; '))) %>% 
    dplyr::left_join(sessionResults$bestRuns.params,.) %>% 
    dplyr::distinct(bid,.keep_all=T) %>% 
    dplyr::select(-c(variable,value))
  
  # find based int method per pol fid
  pivotGrps <- 
    c('testErr.inc.r','compareDiff.inc.r')
  # check if smper present, via testErr.ex.r mean
  testErr.ex.r.mean <- sessionResults$bestRuns %>% 
    dplyr::filter(error_var == 'testErr.ex.r') %>% 
    summarise(mean = mean(error_value))
  
  if (testErr.ex.r.mean != 0) smpper = T else smpper = F
  if (smpper) pivotGrps <- c(pivotGrps,'testErr.ex.r')
  
  bestInts <- sessionResults$bestRuns %>% 
    left_join(bestRuns.params.sum[,c('bid','params_conc')], by='bid') %>% 
    dplyr::group_by(intpol_fid, error_var) %>% 
    slice_min(order_by = error_value)
  
  if (!smpper) { bestInts <- 
    bestInts %>% dplyr::filter(error_var != 'testErr.ex.r') %>% 
    dplyr::filter(error_var %in% c('testErr.inc.r','compareDiff.inc.r')) 
  } else {
    bestInts <- bestInts %>% 
      dplyr::filter(error_var %in% 
                      c('testErr.inc.r','compareDiff.inc.r','testErr.ex.r')) 
  }
  bestIntsFinal <- bestInts %>% 
    dplyr::select(-c(bid,qv,quantile,run_no)) %>%
    tidyr::pivot_wider(names_from = error_var,
                       values_from = c(int_method,params_conc,error_value)) %>% 
    dplyr::mutate_if(is.character,~str_replace_all(.,' with NA = NA','')) %>% 
    dplyr::rename('Site ID' = 'intpol_fid',
                  'Interpolated surface and earlier survey (polygon); method and parameters' = "params_conc_compareDiff.inc.r",
                  'Interpolated surface and later survey (polygon); method and parameters' = "params_conc_testErr.inc.r",
                  'Interpolated surface and earlier survey (polygon): RMSE' = "error_value_compareDiff.inc.r",
                  'Interpolated surface and later survey (polygon): RMSE' = "error_value_testErr.inc.r") %>% 
    { if(smpper) dplyr::rename(.,'Interpolated surface and later survey (noise); method and parameters' = "params_conc_testErr.ex.r",
                            'Interpolated surface and later survey (noise); RMSE' = "error_value_testErr.ex.r") else .} %>% 
    left_join(abDiffs) %>% 
    left_join(terrainA)
  
  # plotting
  
  plots <- list()
  
  # shows relationship between int method, testErr 
  # and terrain metrics
  bestIntsFinal.m <- bestIntsFinal %>% 
    reshape2::melt(id.vars='int_method_testErr.inc.r',
                   measure.vars=
                     setdiff(names(terrainA),'Site ID'))
  
  plots$int_method_testErr.inc.r_and_terrain <- 
    ggplot(bestIntsFinal.m) + 
    geom_boxplot(aes(x = int_method_testErr.inc.r,
                     y = value)) + 
    facet_wrap(~variable,scales='free_y')
  

  # bestInts <- bestInts %>% 
  #   dplyr::filter(`Site ID` %!in% c(21,13,10))
  
  # shows relationship between test RMSE and compare diff rmse
  g1 <- ggplot(bestIntsFinal,
         aes(x = `Interpolated surface and later survey (polygon): RMSE`,
             y = `Interpolated surface and earlier survey (polygon): RMSE`)) + 
    geom_point() +
    stat_smooth(method = 'lm')
  
  g2 <- ggplot(bestIntsFinal,
         aes(x = `Earlier survey and later survey (polygon): RMSE`,
             y = `Interpolated surface and earlier survey (polygon): RMSE`)) + 
    geom_point() +
    stat_smooth(method = 'lm')
  
  rmseScatters <- grid.arrange(g1,g2,ncol=2)
  
  plots$rmseScatters <- rmseScatters
  
  ggplot(bestInts,
         aes(y = `Interpolated surface and earlier survey (polygon): RMSE`,
             x = `Interpolated surface and later survey (noise); RMSE`)) + 
    geom_point() +
    stat_smooth(method = 'lm') + 
    ggrepel::geom_label_repel(aes(label = `Site ID`),
                              nudge_y = 0.01)
  l <- list(plots = plots,
       bestInts = bestInts,
       bestIntsFinal = bestIntsFinal,
       bestInt.l.r = bestInt.l.r)
  return(l)
}

# tagAll <- sessionTag
# tagNew <- paste0(sessionTag,'_bicubic')
bindCubic <- function(tagAll,            # the session tag
                      tagNew,            # the session tag with bicubic appended?
                      outputTag = 'all', # what to append to the bound data output
                      outputDir="/media/mal/working_files/quarry") {
  # name correlations
  nCor <- data.frame(long_name = c("Nearest Neighbor",
                                   "Inverse Distance Weighted",
                                   "Ordinary Kriging",
                                   "GRASS Regularized Splines Tension",
                                   "GRASS Bicubic Spline",
                                   "GRASS Resampled Filter"),
                     short_name = c('nn','idw','ok','gspline',
                                    'gbicubic','gfilter'))
  # find all outputs relating to session tag / 'tagAll'
  f <- list.files(outputDir,pattern=tagAll,
                  full.names = T)
  # removes any outputs with bicubic in filename
  f <- f[!str_detect(f,'bicubic')]
  # finds bicubic outputs relating to the session tag / 'tagNew'
  f.new <- list.files(outputDir,pattern=tagNew,
                      full.names = T)
  
  # copy sources files to backup directory
  print('copying source files to /backup directory...')
  file.copy(f.new,paste0('/media/mal/working_files/quarry/backup/',basename(f.new)))
  file.copy(f,paste0('/media/mal/working_files/quarry/backup/',basename(f)))
  print('done!')
  
  print('binding new data to all other interpolations....')
  # plan("multisession", workers = 6)
  if (length(f) > 0 && length(f.new) > 0 && 
      length(f) == length(f.new) ) {
    # x <- 1
    1:length(f) %>% 
      furrr::future_map(.f = function(x) {
        
        # check pol IDs are the same
        polfid.new <- str_replace(basename(f.new[x]),paste0(tagAll,'_bicubic_'),'')
        polfid <- str_replace(basename(f[x]),paste0(tagAll,'_'),'')
        if (polfid != polfid.new) return(NULL)
        
        print(f.new[x])
        load(f.new[x])
        datNew <- dat
        load(f[x]) # some data saved as 'o' for some reason?
        if (exists('o')) {
          dat <- o
          rm(o)
        }
        
        dat$orig.maps$ras$`GRASS Bicubic Spline` <- 
          datNew$orig.maps$ras$`GRASS Bicubic Spline`
        dat$rmses <- rbind(datNew$rmses,
                           dat$rmses)
        dat$orig.maps$params.m <- 
          rbind(dat$orig.maps$params.m,
                datNew$orig.maps$params.m)
        
        f.out <- str_replace(f[x],tagAll,paste0(tagAll,'_',outputTag))
        save(dat, file=f.out)
        # file.remove(c(f[x],f.new[x]))
        return()
      })
  }
  print('done!')
  print('removing source files from working directory...')
  file.remove(f)
  file.remove(f.new)
  
  m <- str_replace(tagAll,tagAll,paste0(tagAll,'_',outputTag))
  return(m)
}




