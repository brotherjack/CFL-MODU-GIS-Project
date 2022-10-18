#!/usr/bin/env python
from typing import List, Tuple, Optional
from ebird.api import get_species_observations, get_regions
from geojson import Feature, Point, FeatureCollection
from geojson import dump as geojson_dump
from geojson import load as geojson_load
import os


class EbirdManager:
    def __init__(self, api_key, fname=None, targets=[]):
        self.file_name = fname
        self.features = {}
        self.targets = targets
        self.api_key = api_key
        self.geo_json = FeatureCollection(features=[])
        self._raw_pulls = []

    def _is_feature_collection(self, gjson) -> bool:
        if 'type' in gjson.keys() and 'features' in gjson.keys():
            return gjson['type'] == 'FeatureCollection' and hasattr(gjson['features'], 'append')

    def indiv(self, subid):
        if subid in self.features.keys():
            return self.features[subid]['properties']['individuals']
        else:
            return -1

    def checklists_in_geojson(self):
        return [entry['properties']['ebird_subId'] for entry in self.geo_json['features']]

    def set_individuals(self, subid, new_count):
        self.features[subid]['properties']['individuals'] = new_count

    def dedupe(self) -> List[dict]:
        """
        Deduplicates geojson feature list by checking the ebird geojson subid.
        """
        dups = {}
        deduped = []
        for i, feat in enumerate(self.geo_json['features']):
            subid = feat['properties']['ebird_subId']
            if subid in dups.keys():
                if feat == self.geo_json['features'][dups[subid][0]]:
                    dups[subid].append(i)
            else:
                dups[subid] = [i]
                deduped.append(feat)
        return deduped     

    def _setup_features(self):
        for feature in self.geo_json.get('features'):
            subid = feature['properties']['ebird_subId']
            assert(self.indiv(subid) == -1)
            self.features[subid] = feature

    def load(self, fname=None):
        if not fname:
            fname = self.file_name
        
        with open(fname, 'r') as f:
            gjson = geojson_load(f)
            if self._is_feature_collection(gjson):
                self.geo_json = gjson
                self._setup_features()
                print(f"Loaded {len(self.features)} ebird checklist IDs")
            else:
                raise IOError(f"File '{fname}' is not a feature collection!")

    def save(self, fname=None):
        if not fname:
            fname = self.file_name
        
        if os.path.exists(fname):
            with open(fname, 'w') as f:
                geojson_dump(self.geo_json, f)

    def _pull_species_observations(self, region_code, species_codes):
        self._raw_pulls = [] # Reset, if need be
        observations = {}
        
        for species_code in species_codes:
            obs = get_species_observations(self.api_key, species_code, region_code) 
            self._raw_pulls.extend(obs)
            cross_count = 0 # Lumped species in multiple checklists 
            # (eg. "Mottled Duck (Florida)" and "Mallard/Mottled Duck")
            for entry in obs:
                if entry['subId'] in observations.keys():
                    cross_count += 1
                    observations[entry['subId']]['howMany'] += entry['howMany']
                else:
                    observations[entry['subId']] = entry 
            print(f"Pulled {len(obs)} from ebird for {species_code} in {region_code}")
            if cross_count > 0:
                print(f"Found {cross_count} 'lumped' species")

        print(f"Pulled {len(self._raw_pulls)} from ebird for {species_codes} in {region_code}")
        return observations

    def _create_feature_from_observation(self, obs):
        return Feature(
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

    def pull_new_entries(self, region_code, species_codes=None, save=True):
        species_codes = species_codes if species_codes else self.targets
        if not species_codes:
            raise AttributeError(
                "list of ebird 'species_codes' must be passed to method if targets are not set"
            )

        if os.path.exists(self.file_name):
            self.load()
        
        species_obs = self._pull_species_observations(region_code, species_codes)
        gjentries = self.checklists_in_geojson()

        newCount = 0
        for subid, obs in species_obs.items():
            subid_count = self.indiv(subid)
            if subid not in gjentries:
                newCount += 1
                self.geo_json.features.append(self._create_feature_from_observation(obs))
                self.features[subid] = self.geo_json.features[-1]
            else:
                if subid_count != obs['howMany']:
                    print(
                        f"Count for {subid} has changed", 
                        f"from {obs['howMany']} to {subid_count}"
                    )
                    self.set_individuals(subid, subid_count)
                    index = gjentries.index(subid)
                    assert(self.geo_json.features[index]['properties']['ebird_subId'] == subid)
                    self.geo_json.features[index] = self._create_feature_from_observation(obs)
        
        if save:
            self.save()

        return newCount


if __name__ == '__main__':
    API_KEY = os.environ.get("API_KEY")

    if not API_KEY:
        raise IOError("Environmental variable 'API_KEY' must be set")

    region_code = 'US-FL-095'
    species_codes = ['motduc', 'x00422', 'motduc1', 'y00632']

    ebird_man = EbirdManager(API_KEY, "ebird.geojson", targets=species_codes)
    new_count = ebird_man.pull_new_entries(region_code, save=True)
    print(f"Saved {new_count} new entries")

