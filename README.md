# Local variation in COVID excess mortality across Italy

This is a project being developed by Nat Henry, nathaniel.henry@ndm.ox.ac.uk

## Downloading Data

Mortality data, population data, and census-based covariates were downloaded 
from IStat, the Italian Statistical Authority.

#### Mortality

IStat has periodically released detailed tabulations of all-cause mortality by
municipality to help explore the effects of COVID on mortality. I accessed the
[COVID mortality web page](https://www.istat.it/it/archivio/240401) on August 18,
2020 and downloaded the "Dataset con i decessi giornalieri" (Dataset with daily 
deaths) which had last been updated on August 10, 2020. Direct link to the 
mortality dataset as of August 18, 2020: 
https://www.istat.it/it/files//2020/03/Dataset-decessi-comunali-giornalieri-e-tracciato-record_al30giugno.zip

#### Population

Population data was downloaded from the [IStat Data Portal](https://dati.istat.it) 
under the header "Population and Households > Population > Resident Population 
on 1st January > All municipalities". I created a free account for the download 
and interacted with the Italian-language version of the website. I then used the
following settings to download tabulated population data:

- **Layout:**
  - Filters: Demographic variable, sex, year, marriage status
  - Vertical dimensions: Territory (country/province/municipality), Age
  - Horizontal dimensions: (None)

- **Filters:** I downloaded separate datasets for each sex (male/female) and each year
  (2015/2016/2017/2018/2019/2020) available in the COVID excess mortality data.

- **Export:** I exported to a CSV, using "custom format" export options:
  - Included both codes and labels by field
  - English language
  - Comma separated

#### Covariates

Covariates were downloaded by commune for all available years since 2015:
- Total fertility rate (2015-2018), by province
- Proportion of targeted families reached by various social services (2015-2017), by province
- Q1 unemployment by sex, all age groups (2015-2020), by province


For more recent years where a covariate was missing, the average of the last two
available years of data was taken

#### Shapefile

A commune-level shapefile was downloaded from the [IStat website](https://www.istat.it/it/archivio/222527).
Direct link to 2020 shapefile as of August 18, 2020: 
- Most detailed (for analysis): http://www.istat.it/storage/cartografia/confini_amministrativi/non_generalizzati/Limiti01012020.zip
- Less detailed (for mapping): http://www.istat.it/storage/cartografia/confini_amministrativi/generalizzati/Limiti01012020_g.zip


## License

This repository operates under the GNU General Public License version 3.0. For 
more details, see the `LICENSE` file in this repository.
