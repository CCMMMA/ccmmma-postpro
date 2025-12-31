Run the following script to have data playgroud.

```bash
cd data
prods="wrf5 ww33 rms3 wcm3 aiq3"
mkdir history archive
for prod in $prods
do
  echo $prod
  wget https://data.meteo.uniparthenope.it/files/$prod/d03/archive/2026/01/01/${prod}_d03_20260101Z1200.nc -O archive/${prod}_d03_20260101Z1200.nc
  wget https://data.meteo.uniparthenope.it/files/$prod/d03/history/2026/01/01/${prod}_d03_20260101Z1200.nc -O history/${prod}_d03_20260101Z1200.nc
done
```

