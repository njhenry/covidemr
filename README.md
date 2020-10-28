# Local variation in COVID excess mortality across Italy

This is a project being developed by Nat Henry, nathaniel.henry@ndm.ox.ac.uk

## Downloading Data

Mortality data, population data, and census-based covariates were downloaded
from IStat, the Italian Statistical Authority.

### Mortality

IStat has periodically released detailed tabulations of all-cause mortality by
municipality to help explore the effects of COVID on mortality. I accessed the
[COVID mortality web page](https://www.istat.it/it/archivio/240401) on August 18,
2020 and downloaded the "Dataset con i decessi giornalieri" (Dataset with daily
deaths) which had last been updated on August 10, 2020. Direct link to the
mortality dataset as of August 18, 2020, updated through June 30:
https://www.istat.it/it/files//2020/03/Dataset-decessi-comunali-giornalieri-e-tracciato-record_al30giugno.zip

### Population

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

### Covariates:

#### IStat Covariates

Covariates were downloaded by province for all available years since 2015:
- Total fertility rate (2015-2018), by province
- Proportion of targeted families receiving at-home social services (2015-2017),
  by province
- Annual unemployment by sex, ages 15 and above (2015-2019), by province
- Proportion of households with taxable income under 10k Euros (2015-2018), by
  commune (aggregated to province)
- Average taxable income across all households (2015-2018), by commune
  (aggregated to province)

#### Elevation

I downloaded a digital elevation map at 15 arc-second resolution from the
US Geological Survey Earth Explorer portal:

1. Create a free web account with the [USGS EROS portal](https://ers.cr.usgs.gov/)
2. Navigate to the USGS EarthExplorer: https://earthexplorer.usgs.gov/
3. Create a bounding box containing Italy: I used the coordinates bounded by
   [6.5 deg E, 35 deg N] and [20 deg E, 48 deg N]
4. From the "Data Sets" tab, select the "Digital Elevation > GMTED2010" product.
   This searches only for the [USGS Global Multi-Resolution Terrain Elevation
   Data 2010](https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-global-multi-resolution-terrain-elevation)
   product.
5. Download and unzip TIFF file.

#### Healthcare access

I downloaded a gridded dataset showing travel time to the nearest health
facility by motor vehicle (see [Weiss et al 2020](https://www.nature.com/articles/s41591-020-1059-1)) using a direct download from the Malaria Atlas Project [Data
Explorer](https://malariaatlas.org/explorer/#/).

#### Weekly temperature

I downloaded point estimates of daily temperature from the three most populous
pixels in each province using the [Meteostat API](https://dev.meteostat.net/docs/api/).
This required signing up for a free API key
Because temperature data was not available for all time periods and locations, I
used the following process to fill temperatures:

1. Find average weekly temperature across all observed days by year, week, and
   observed location.
2. Interpolate by week: in cases where 3 weeks or fewer were missing between
   observed data, interpolate temperature from neighboring weeks in a given
   observation location.
3. For each province, year, and week, average available observations across the
   three sampled sites
4. In cases where all observations for a province were missing for a given
   province-year-week, fill using temperature observations from neighboring
   provinces with a similar [elevation](http://tinitaly.pi.ingv.it/) and
   [level of solar exposure](https://solargis.com/maps-and-gis-data/download/italy).


### Shapefile

A commune-level shapefile was downloaded from the [IStat website](https://www.istat.it/it/archivio/222527).
Direct link to 2020 shapefile as of August 18, 2020:
- Most detailed (for analysis): http://www.istat.it/storage/cartografia/confini_amministrativi/non_generalizzati/Limiti01012020.zip
- Less detailed (for mapping): http://www.istat.it/storage/cartografia/confini_amministrativi/generalizzati/Limiti01012020_g.zip


## License

This repository operates under the GNU General Public License version 3.0. For
more details, see the `LICENSE` file in this repository.
