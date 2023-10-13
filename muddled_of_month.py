mall_df = gpd.read_file('mall.geojson')
d = (dt.datetime.today() - dt.timedelta(weeks=4)).date()
mall_df = mall_df[mall_df.observation_date >= f"{d.year}-{d.month}-{d.day}"]
mall_df = mall_df.assign(species="mallard")


df = gpd.read_file('ebird.geojson')
df = df[df.observation_date >= f"{d.year}-{d.month}-{d.day}"]
df = df.assign(species="mottled duck")

df = pd.concat([df, mall_df])
t = dt.datetime.today().date()
df.to_file(f"/home/thellinger/Downloads/modu_x_mall_{d.year}_{d.month:02}_{d.day:02}__{t.year}_{t.month:02}_{t.day:02}.geojson", driver='GeoJSON')
