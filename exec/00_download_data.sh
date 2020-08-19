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

echo "=========================================" &&
echo "*** Data will be saved to $1 ***" &&
cd $1 && mkdir shp && mkdir deaths && cd shp && 
wget http://www.istat.it/storage/cartografia/confini_amministrativi/non_generalizzati/Limiti01012020.zip &&
unzip Limiti01012020.zip && 
wget http://www.istat.it/storage/cartografia/confini_amministrativi/generalizzati/Limiti01012020_g.zip &&
unzip Limiti01012020_g.zip &&
echo "*** - Detailed and generalized shapefiles saved to $1/shp ***" && 
cd ../deaths && 
wget https://www.istat.it/it/files//2020/03/Dataset-decessi-comunali-giornalieri-e-tracciato-record_al30giugno.zip
unzip Dataset-decessi-comunali-giornalieri-e-tracciato-record_al30giugno.zip &&
echo "*** - Deaths saved to $1/deaths ***" &&
echo "====================DONE====================";