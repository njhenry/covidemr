## -----------------------------------------------------------------------------
##
## 00: Italy data prep: Download data
##
## Download raw shapefile and death data from Istat
## For more information, see README at https://github.com/njhenry/covidemr/
##
## Execution: Takes one positional argument, location of raw data folder
## Example:
##   ./00_download_data.sh /path/to/input/data/folder/
##
## -----------------------------------------------------------------------------

echo "======================================================================" &&
echo "" &&
echo "*** Data will be saved to $1 ***" &&
cd $1 && mkdir -p shp deaths pop_raster covars/raster covid_deaths &&
cd shp &&
wget http://www.istat.it/storage/cartografia/confini_amministrativi/non_generalizzati/Limiti01012020.zip &&
wget http://www.istat.it/storage/cartografia/confini_amministrativi/generalizzati/Limiti01012020_g.zip &&
wget https://www.istat.it/storage/codici-unita-amministrative/Elenco-codici-statistici-e-denominazioni-delle-unita-territoriali.zip\
  -O location_codes.zip &&
unzip Limiti01012020.zip && unzip Limiti01012020_g.zip && unzip location_codes.zip &&
cp Elenco*/*.csv ./location_codes.csv &&
rm -r Elenco*/ && rm location_codes.zip && rm Limiti01012020.zip && rm Limiti01012020_g.zip &&
echo "***   - Detailed and generalized shapefiles saved to $1/shp ***" &&
cd ../deaths &&
wget https://www.istat.it/it/files//2020/03/Dataset-decessi-comunali-giornalieri-e-tracciato-record_22ottobre2020.zip\
  -O daily_deaths.zip &&
unzip daily_deaths.zip &&
echo "***   - Deaths saved to $1/deaths ***" &&
echo "*** Downloading population raster: ***" &&
cd ../pop_raster &&
WP_BASE=ftp://ftp.worldpop.org.uk/GIS/Population/Global_2000_2020_1km &&
curl $WP_BASE/2015/ITA/ita_ppp_2015_1km_Aggregated.tif -o ita_2015.tif &&
curl $WP_BASE/2015/ITA/ita_ppp_2015_1km_Aggregated.tif -o ita_2015.tif &&
curl $WP_BASE/2016/ITA/ita_ppp_2016_1km_Aggregated.tif -o ita_2016.tif &&
curl $WP_BASE/2017/ITA/ita_ppp_2017_1km_Aggregated.tif -o ita_2017.tif &&
curl $WP_BASE/2018/ITA/ita_ppp_2018_1km_Aggregated.tif -o ita_2018.tif &&
curl $WP_BASE/2019/ITA/ita_ppp_2019_1km_Aggregated.tif -o ita_2019.tif &&
curl $WP_BASE/2020/ITA/ita_ppp_2020_1km_Aggregated.tif -o ita_2020.tif &&
echo "*** Downloading COVID deaths from healthdata.org ***"
cd ../covid_deaths &&
wget https://ihmecovid19storage.blob.core.windows.net/archive/2020-10-15/ihme-covid19.zip &&
unzip ihme-covid19.zip &&
mv 2020_*/* ./ && rm -r 2020_* && rm ihme-covid19.zip &&
echo "*** Copying Meteostat API key to date-specific folder ***" &&
cd ../covars/ &&
cp ../../meteostat_api_key.txt ./ &&
mkdir meteostat_cache_dir &&
echo "*** Downloading raster covariates: ***" &&
echo "***   - Access to hospitals ***" &&
cd raster &&
curl "https://malariaatlas.org/geoserver/ows?service=CSW&version=2.0.1&request=DirectDownload&ResourceId=Explorer:2020_motorized_travel_time_to_healthcare" \
  -o healthcare_access.zip &&
unzip healthcare_access.zip &&
rm healthcare_access.zip &&
mv 2020_motorized_travel_time_to_healthcare.geotiff healthcare_access.tif &&
echo "" &&
echo "================================ DONE ================================";
