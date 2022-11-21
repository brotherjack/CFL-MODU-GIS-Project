#!/usr/bin/env python
from typing import List, Tuple, Optional
import argparse
from collections import OrderedDict
from IPython.terminal.embed import InteractiveShellEmbed
from ebird.api import get_species_observations, get_regions
from geojson import Feature, Point, FeatureCollection
from geojson import dump as geojson_dump
from geojson import load as geojson_load
import numpy as np
import geopandas as gpd
import pandas as pd
import logging, os, re, unicodedata, uuid


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


CHECK_MARK = unicodedata.lookup('WHITE HEAVY CHECK MARK')
NOPE_MARK =  unicodedata.lookup('NO ENTRY SIGN')
THUMBS_UP = u'\U0001f44d'
SHRUGGIE = "¯\_(ツ)_/¯"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if logger.hasHandlers():
    logger.handlers[0].setFormatter(logging.Formatter(logging.BASIC_FORMAT))
else:
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
    logger.addHandler(handler)

try:
    get_ipython
except NameError:
    banner=exit_msg=''
banner = """
                   ________________    __  _______  ____  __  __
                  / ____/ ___/ ___/   /  |/  / __ \/ __ \/ / / /
                 / /    \__ \ \__\   / /|_/ / / / / / / / / / /
                / /___ ___/ /__/ /  / /  / / /_/ / /_/ / /_/ /
                \____//____/____/  /_/  /_/\____/_____/\____/

                         +-+-+ +-+-+ +-+-+-+-+-+-+-+-+
                         |I|T| |I|S| |D|U|C|K|T|I|M|E|
                         +-+-+ +-+-+ +-+-+-+-+-+-+-+-+

"""
exit_msg = """
            +-+-+-+ +-+-+-+ +-+-+-+-+-+ +-+-+-+-+-+-+-+
            |S|e|e| |Y|o|u| |S|p|a|c|e| |D|u|c|k|b|o|i|
            +-+-+-+ +-+-+-+ +-+-+-+-+-+ +-+-+-+-+-+-+-+
"""
ipshell = InteractiveShellEmbed(
    banner1=f"{bcolors.FAIL} {banner} {bcolors.ENDC}", 
    exit_msg=f"{bcolors.FAIL} {exit_msg} {bcolors.ENDC}"
)


class EbirdManager:
    def __init__(self, api_key, fname=None, targets=[]):
        self.file_name = fname
        self.features = {}
        self.targets = targets
        self.api_key = api_key
        self.geo_json = FeatureCollection(features=[])
        self._ebird_regex = re.compile("S[0-9]{9}") # Regex for ebird checklist IDs
        self._raw_pulls = []
        self._trapping_area_cache = None
        self.modu_site_reports = None
    
    def _modu_count_column_mapper(self, col):
        col = col.lower().replace(" ", "_")
        if "longitude" in col:
            return "longitude"
        elif "latitude" in col:
            return "latitude"
        elif "modu" in col:
            return "modu_count"
        elif "whib" in col:
            return "whib_count"
        elif "mudu" in col:
            return "mudu_count"
        elif "grgo" in col:
            return "grgo_count"
        else:
            return col

    def _clean_modu_site_df_columns(self, df) -> pd.DataFrame:
        df = df.rename(mapper=self._modu_count_column_mapper, axis=1)
        df.modu_count = df.modu_count.convert_dtypes()
        df.whib_count = df.whib_count.convert_dtypes()
        df.mudu_count = df.mudu_count.convert_dtypes()
        df.grgo_count = df.grgo_count.convert_dtypes()
        df.scout_date = df.scout_date.dt.date
        df.added = df.added.astype(bool)
        df.s3 = df.s3.str.slice(0,1).astype(int)
        return df

    def is_checklist_id_valid(self, checklist_id) -> bool:
        """Checks if ebird checklist ID is formated correctly. 
        """
        if type(checklist_id) == str:
            return type(self._ebird_regex.match(checklist_id)) == re.Match 
        else:
            return False

    def _is_feature_collection(self, gjson) -> bool:
        if 'type' in gjson.keys() and 'features' in gjson.keys():
            return gjson['type'] == 'FeatureCollection' and hasattr(gjson['features'], 'append')

    def indiv(self, subid):
        """
        Returns the number of individuals for an ebird checklist, 
        or -1 if the checklist ID is not found
        """
        if subid in self.features.keys():
            return self.features[subid]['properties']['individuals']
        else:
            return -1

    def checklists_in_geojson(self) -> list:
        """Returns checklist ID's in geojson data
        """
        return [entry['properties']['ebird_subId'] for entry in self.geo_json['features']]

    def set_individuals(self, subid, new_count):
        """Changes the count of a checklist
        """
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
        """
        Loads the geojson file at fname (or in self.file_name). Sets up 
        this individual's feature hash. 
        """
        if not fname:
            fname = self.file_name
        
        with open(fname, 'r') as f:
            gjson = geojson_load(f)
            if self._is_feature_collection(gjson):
                self.geo_json = gjson
                self._setup_features()
                logger.info(f"Loaded {len(self.features)} ebird checklist IDs")
            else:
                raise IOError(f"File '{fname}' is not a feature collection!")

    def save(self, fname=None) -> None:
        """Saves gejson object into a geojson file at (fname or self.file_name)
        """
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
            logging.info(f"Pulled {len(obs)} from ebird for {species_code} in {region_code}")
            if cross_count > 0:
                logging.debug(f"Found {cross_count} 'lumped' species")

        logging.debug(f"Pulled {len(self._raw_pulls)} from ebird for {species_codes} in {region_code}")
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

    def pull_new_entries(self, region_code, species_codes=None, save=True) -> int:
        """
        Takes an ebird `region_code` (eg. "US-FL-095") and a list of ebird
        `species_codes` (eg. ['motduc', 'x00422', 'motduc1', 'y00632']),
        defaults to `self.targets`. This method then uses this data to make
        an API call to ebird to download the latest observations for the 
        `species_codes`, which are treated as equivalent (ie. "mottled duck" 
        and "hybrid mottled x mallard duck") are counted together. New entries
        are appended to `self.geo_json` and this object and `self.features` 
        are upated as they are in the ebird database. If `save` is true,
        `self.save()` is called.

        Returns the number of new entries.
        """
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
                    logging.warning(
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

    def import_modu_site_reports(self, fname):
        """Imports modu site reports from Google forms XLSX
        """
        self.modu_site_reports = self._clean_modu_site_df_columns(
            pd.read_excel(fname)
        )

    def import_survey_sites(self, fname="survey_sites.gpkg", layer='orlando_parks'):
        """Imports survey sites from geopackage
        """
        self.survey_sites = gpd.read_file(fname, layer=layer)
        for introw in ['MODU_COUNT', 'MUDU_COUNT', 'WHIB_COUNT', 'GRGO_COUNT', 'AREA', 'OBJECTID']:
            setattr(self.survey_sites, introw, getattr(self.survey_sites, introw).convert_dtypes())
        self.survey_sites.replace({np.nan: None}, inplace=True)

    def export_survey_sites(self, fname="survey_sites.gpkg", layer='orlando_parks'):
        self.survey_sites.to_file(fname, layer=layer, driver="GPKG")
        
    def names_of_duplicated_sites(self):
        return sorted(list(self.survey_sites[self.survey_sites.duplicated()].NAME))

    def import_scouting_areas(self, fname="scouting_regions.shp"):
        self.scouting_areas = gpd.read_file(fname)

    def _global_id_is_unique(self, entry):
        if len(entry) != 1:
            # TODO: Should this be a critical error?
            raise IndexError(
                    f"There should only be 1 entry with global ID "
                    f"'{entry.GlobalID.iloc[0]}' however there are "
                    f"{len(entry)} entries with that ID"
                )
        else:
            return True

    def find_scouting_area_for_site(self, global_id):
        entry = self.survey_sites[self.survey_sites.GlobalID == global_id]

        try:
            if self._global_id_is_unique(entry):
                entry = entry.iloc[0]
            else:
                return False
        except IndexError as ie:
            logger.error(ie)

        for ind in self.scouting_areas.index:
            if entry.geometry.within(self.scouting_areas.loc[ind, 'geometry']):
                return self.scouting_areas.loc[ind, "id"]

        return None

    def validation_pass(self, msg, level="info", suppress=False):
        if not suppress:
            getattr(logger, level)(f"{bcolors.OKGREEN} {msg} {bcolors.ENDC}")

    def validation_fail(self, msg, level="warning", suppress=False):
        if not suppress:
            getattr(logger, level)(f"{bcolors.FAIL} {msg} {bcolors.ENDC}")

    def verify_global_ids_unique(self):
        valid = True
        duplicated_global_ids = set(
            self.survey_sites[self.survey_sites.GlobalID.duplicated()].GlobalID
        )

        if duplicated_global_ids:
            valid = False
            for duplicated_global_id in duplicated_global_ids:
                ids = self.survey_sites[self.survey_sites.GlobalID == duplicated_global_id]
                self.validation_fail(
                    f"{duplicated_global_id} is the Global ID for {list(ids.index)}"
                )

        if valid:
            self.validation_pass(
                f"{bcolors.OKGREEN} {CHECK_MARK} Global IDs are unique "
                f"{THUMBS_UP} {bcolors.ENDC}"
            )
        return valid

    def verify_trapping_area(self, ind, supress_log=False):
        row = self.survey_sites.loc[ind, :]
        found_area = self.find_scouting_area_for_site(row.GlobalID)
        if row.AREA:
            if str(row.AREA) != str(found_area):
                self.validation_fail(
                    f"Trapping area for {row.NAME} - {row.GlobalID} is "
                    f"incorrect. Is marked as {row.AREA} should be "
                    f"{found_area}",
                    suppress=supress_log
                )
                self._trapping_area_cache = str(found_area)
                return False
            else:
                self.validation_pass(
                    f"Trapping area for {row.NAME} is correct as {row.AREA}",
                    level="debug",
                    suppress=supress_log
                )
        else:
            if found_area:
                self.validation_fail(
                    f"Trapping area for {row.NAME} is blank but should "
                    f"be {found_area}",
                    suppress=supress_log
                )
                self._trapping_area_cache = str(found_area)
                return False
            else:
                self.validation_pass(
                    f"Trapping area for {row.NAME} is blank "
                    f"{bcolors.UNDERLINE}and should be so",
                    level="debug",
                    suppress=supress_log
                )
        return True

    def verify_trapping_areas(self):
        valid = True
        for ind in self.survey_sites.index:
            if not self.verify_trapping_area(ind):
                valid = False

        if valid:
            self.validation_pass(
                f"{CHECK_MARK} All survey sites have correct trapping "
                f"areas {THUMBS_UP}"
            )
        else:
            self.validation_fail(
                f"{NOPE_MARK} At least some survey sites are wrong"
            )

        return valid

    def verify_survey_sites(self, only=[]):
        VERIFICATIONS = set([
            "GLOBAL_IDS_UNIQUE",
            "TRAPPING_AREAS"
        ]).difference(set([s.upper() for s in only]))

        return all([
            getattr(self, f"verify_{verify_suffix.lower()}")()
            for verify_suffix in VERIFICATIONS
        ])

    def add_global_id(self, ind, overwrite=False):
        lname = {self.survey_sites.loc[ind, 'NAME']}
        if self.survey_sites.loc[ind].notnull()['GlobalID']:
            old_id = self.survey_sites.loc[ind, 'GlobalID']
            if overwrite:
                new_id = "{" + f"{uuid.uuid4()}".upper() +"}"
                logger.info(
                    f"Overwritting GlobalID for {lname} "
                    f"from {old_id} to {new_id}"
                )
            else:
                logger.warning(
                    f"Will not update {lname} GlobalID as it already exsists"
                    f"as {old_id} and overwrite is set as 'False'"
                )
        else:
            new_id = "{" + f"{uuid.uuid4()}".upper() +"}"
            self.survey_sites.loc[ind, 'GlobalID'] = new_id
            logger.info(f"Added GlobalID {new_id} to {lname}")

    def add_global_ids(self):
        for ind in self.survey_sites.index:
            if self.survey_sites.loc[ind].isnull()['GlobalID']:
                self.add_global_id(ind)

    def correct_trapping_areas(self, supress_log=True):
        for ind in self.survey_sites.index:
            if not self.verify_trapping_area(ind, supress_log):
                name, area, gid = self.survey_sites.loc[ind, ["NAME", "AREA", "GlobalID"]]
                self.survey_sites.loc[ind, "AREA"] = self._trapping_area_cache
                logger.info(
                    f"Corrected trapping area for {name} - {gid}, "
                    f"from '{area}' to '{self._trapping_area_cache}'"
                )
        logger.info(f"Corrections completed {THUMBS_UP}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--level',
        '-l',
        required=False,
        type=str.upper,
        choices=["DEBUG", "INFO", "WARNING", "WARN", "CRITICAL", "ERROR"]
    )

    parser.add_argument('--pull', '-p', action="store_true", required=False)
    parser.add_argument('--verify', '-v', action="store_true", required=False)
    parser.add_argument('--interpreter', '-i', action="store_true", required=False)
    parser.add_argument('--correct', '-c', action="store_true", required=False)
    args = parser.parse_args()

    if args.level:
        logger.setLevel(getattr(logging, args.level))

    API_KEY = os.environ.get("API_KEY")

    if not API_KEY:
        raise IOError("Environmental variable 'API_KEY' must be set")

    region_code = 'US-FL-095'
    species_codes = ['motduc', 'x00422', 'motduc1', 'y00632']

    ebird_man = EbirdManager(API_KEY, "ebird.geojson", targets=species_codes)

    if args.pull:
        logger.info("Pulling data...")
        new_count = ebird_man.pull_new_entries(region_code, save=True)
        logger.info(f"Saved {new_count} new entries")

    ebird_man.import_survey_sites()
    ebird_man.import_scouting_areas()

    if args.verify:
        ebird_man.verify_survey_sites()

    if args.correct:
        ebird_man.correct_trapping_areas()
        ebird_man.export_survey_sites()

    if args.interpreter:
        ipshell()
