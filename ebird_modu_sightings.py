#!/usr/bin/env python
from ebird.api import get_species_observations, get_regions
from geojson import Feature, Point, FeatureCollection
from geojson import dump as geojson_dump
from geojson import load as geojson_load
import os
import sys
import datetime as dt


class EbirdManager:
    def __init__(self, api_key, fname=None):
        self.file_name = fname
        self.sub_ids = []
        self.api_key = api_key
        self.geo_json = FeatureCollection(features=[])

    def _is_feature_collection(self, gjson):
        if 'type' in gjson.keys() and 'features' in gjson.keys():
            return gjson['type'] == 'FeatureCollection' and hasattr(gjson['features'], 'append')

    def load(self, fname=None):
        if not fname:
            fname = self.file_name
        
        with open(fname, 'r') as f:
            gjson = geojson_load(f)
            if self._is_feature_collection(gjson):
                return gjson
            else:
                raise IOError(f"File '{fname}' is not a feature collection!")

    def save(self, fname=None):
        if not fname:
            fname = self.file_name
        
        if os.path.exists(fname):
            with open(fname, 'w') as f:
                geojson_dump(self.geo_json, f)

    def pull_new_entries(self, species_codes, region_code, save=True):
        if os.path.exists(self.file_name):
            self.geo_json = self.load()
        
            for feature in self.geo_json.get('features'):
                self.sub_ids.append(feature['properties']['ebird_subId'])
            print(f"Loaded {len(self.sub_ids)} ebird checklist IDs")
        
        modu_obs = get_species_observations(self.api_key, species_codes, region_code)
        print(f"Pulled {len(modu_obs)} from ebird for {species_codes} in {region_code}")
        
        newCount = 0
        for obs in modu_obs:
            if obs["subId"] in self.sub_ids:
                continue
            
            newCount += 1
            self.geo_json.features.append(
                Feature(
                    geometry=Point((obs['lng'], obs['lat'])), properties={
                    "ebird_locid": obs["locId"],
                    "observation_date": obs["obsDt"],
                    "individuals": obs["howMany"],
                    "ebird_valid": obs["obsValid"],
                    "ebird_reviewed": obs["obsReviewed"],
                    "locationPrivate": obs["locationPrivate"],
                    "ebird_subId": obs["subId"],
                    "ebird_location_name": obs["locName"]
                })
            )
    
        with open(self.file_name, 'w') as f:
            geojson_dump(self.geo_json, f)
            print(f"Saved {newCount} new entries")


if __name__ == '__main__':
    API_KEY = os.environ.get("API_KEY")
    
    if not API_KEY:
        raise IOError("Environmental variable 'API_KEY' must be set")

    region_code = 'US-FL-095'
    species_code = 'motduc'
    
    ebird_man = EbirdManager(API_KEY, "ebird.geojson")
    ebird_man.pull_new_entries(species_code, region_code)

