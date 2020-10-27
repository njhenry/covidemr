## -----------------------------------------------------------------------------
##
## 00: Italy data prep: Download data
##
## Download raw shapefile and death data from Istat
## For more information, see README at https://github.com/njhenry/covidemr/
##
## Execution: Takes one positional argument, location of raw data folder
## Example:
##   ./00_download_data.sh /path/to/input/folder/
##
## -----------------------------------------------------------------------------

echo "======================================================================" &&
echo "" &&
echo "*** Data will be saved to $1 ***" &&
cd $1 && mkdir shp && mkdir deaths && mkdir pop_raster && mkdir covars/raw && cd shp &&
wget http://www.istat.it/storage/cartografia/confini_amministrativi/non_generalizzati/Limiti01012020.zip &&
unzip Limiti01012020.zip &&
wget http://www.istat.it/storage/cartografia/confini_amministrativi/generalizzati/Limiti01012020_g.zip &&
unzip Limiti01012020_g.zip &&
echo "***   - Detailed and generalized shapefiles saved to $1/shp ***" &&
cd ../deaths &&
wget https://www.istat.it/it/files//2020/03/Dataset-decessi-comunali-giornalieri-e-tracciato-record_al30giugno.zip
unzip Dataset-decessi-comunali-giornalieri-e-tracciato-record_al30giugno.zip &&
echo "***   - Deaths saved to $1/deaths ***" &&
echo "*** Downloading population raster: ***" &&
cd ../pop_raster &&
WP_BASE=ftp://ftp.worldpop.org.uk/GIS/Population/Global_2000_2020_1km &&
wget $WP_BASE/2015/ITA/ita_ppp_2015_1km_Aggregated.tif -O ita_2015.tif --no-passive &&
wget $WP_BASE/2016/ITA/ita_ppp_2016_1km_Aggregated.tif -O ita_2016.tif --no-passive &&
wget $WP_BASE/2017/ITA/ita_ppp_2017_1km_Aggregated.tif -O ita_2017.tif --no-passive &&
wget $WP_BASE/2018/ITA/ita_ppp_2018_1km_Aggregated.tif -O ita_2018.tif --no-passive &&
wget $WP_BASE/2019/ITA/ita_ppp_2019_1km_Aggregated.tif -O ita_2019.tif --no-passive &&
wget $WP_BASE/2020/ITA/ita_ppp_2020_1km_Aggregated.tif -O ita_2020.tif --no-passive &&
echo "*** Copying Meteostat API key to date-specific folder ***" &&
cd ../covars/ &&
cp ../../meteostat_api_key.txt ./ &&
mkdir meteostat_cache_dir &&
echo "*** Downloading raster covariates: ***" &&
echo "***   - 1) Access to hospitals ***" &&
cd raw &&
wget "https://malariaatlas.org/geoserver/ows?service=CSW&version=2.0.1&request=DirectDownload&ResourceId=Explorer:2020_motorized_travel_time_to_healthcare" \
  -O healthcare_access.zip &&
unzip healthcare_access.zip &&
rm healthcare_access.zip &&
mv 2020_motorized_travel_time_to_healthcare.geotiff healthcare_access.tif &&
echo "***   - 2) Elevation ***" &&
wget http://tinitaly.pi.ingv.it/AJK764GHJ0987NBV/TINITALY_image.zip -O italy_DEM.zip &&
unzip italy_DEM.zip -d italy_DEM &&
rm italy_DEM.zip &&
echo "***   - 3) Average monthly temperature ***" &&
wget https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.04/cruts.2004151855.v4.04/tmp/cru_ts4.04.2011.2019.tmp.dat.gz &&
gzip -d cru_ts4.04.2011.2019.tmp.dat.gz &&
wget https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.04/cruts.2004151855.v4.04/tmp/cru_ts4.04.2011.2019.tmp.dat.nc.gz &&
gzip -d cru_ts4.04.2011.2019.tmp.dat.nc.gz &&
wget https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.04/cruts.2004151855.v4.04/tmp/cru_ts4.04.2011.2019.tmp.stn.gz &&
gzip -d cru_ts4.04.2011.2019.tmp.stn.gz &&
echo "" &&
echo "================================ DONE ================================";
