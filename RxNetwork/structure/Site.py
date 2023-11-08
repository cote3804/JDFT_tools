from __future__ import annotations
from structure.Structure import Structure
from data.Data_Parser import Data_Parser
from helpers.conversions import bias_float_to_str


class Site:
    def __init__(self, sites_dict:dict):
        self.sites_data = sites_dict
    
    def get_site_data(self, site_code):
        return self.sites_data[site_code]

    def __repr__(self):
        print_string = "sites: "
        for site in self.sites_data.keys():
            print_string += f"{site} "
        return print_string

    def get_min_site(self):
        E_min = 0
        min_site = ''
        for site in self.sites_data.keys():
            if "converged" in self.sites_data[site] and self.sites_data[site]["converged"] == True and self.sites_data[site]["final_energy"] < E_min:
                min_site = site
                E_min = self.sites_data[site]["final_energy"]
            else:
                continue
        return min_site

    def unique_clusters(self, intermediate:str):
        # This method returns a list of unique clusters
        # It does this by comparing the current cluster with the other clusters in the list
        # If the current cluster is the same as any of the other clusters, it is not added to the list
        # If the current cluster is not the same as any of the other clusters, it is added to the list
        unique_clusters = []
        for site in self.sites_data.keys():
            if "contcar" in self.sites_data[site].keys(): # if no structure data, ignore site
                cluster = Cluster(self.sites_data, intermediate, site)
            elif "contcar" not in self.sites_data[site].keys():
                continue
            neighbors = cluster.neighbors()
            # this if statement checks if the current site cluster already exists in the unique_clusters list
            if not any([cluster.check_same_cluster(unique_cluster) for unique_cluster in unique_clusters]):
                if len(neighbors) > 0:
                    unique_clusters.append(cluster)
                else:
                    continue
            else:
                continue
        return unique_clusters
    
    @classmethod
    def get_site(cls, data_path, filename, surface, intermediate, bias):
        # instantiatie a Site object from string inputs alone
        bias = bias_float_to_str(bias)
        intermediate = intermediate.strip("*")
        data = Data_Parser(data_path, filename=filename)
        sites_data = data.get_sites_data(surface, bias, intermediate)
        return cls(sites_data)




class Cluster(Site):
    def __init__(self, sites_data:dict, intermediate:str, site:str):
        self.sites_data = sites_data
        self.site_data = sites_data[site]
        self.structure = Structure(self.site_data["contcar"])
        self.intermediate = intermediate

    def __repr__(self): # cluster object prints neighbors of the adsorbate
        neighbors = self.neighbors()
        return f"{neighbors}"
    
    def neighbors(self) -> list:
        neighbors = self.structure.intermediate_neighbors(self.intermediate)
        return neighbors

    def check_in_cluster(self, cluster:Cluster) -> bool:
        # method that compares the current cluster object with another one passed into the method.
        # If both clusters share at least one atom, the cluster is the same.
        neighbors = self.neighbors()
        other_neighbors = cluster.neighbors()
        match_list = [neighbor in neighbors for neighbor in other_neighbors] 
        # this line iterates through other_neighbors and returns a boolean value for each
        # element depending on whether it is in neighbors. This ends as a list of booleans that
        # can be further checked with any() to check overlap of cluster lists
        return any(match_list)
    
    def check_same_cluster(self, cluster:Cluster) -> bool:
        # method that compares the current cluster object with another one passed into the method.
        # If both clusters share at least one atom, the cluster is the same.
        neighbors = set(self.neighbors())
        other_neighbors = set(cluster.neighbors())
        are_same = neighbors == other_neighbors
        return are_same
    
    def ground_site(self):
        sites_in_cluster = []
        for site in self.sites_data.keys():
            if "contcar" not in self.sites_data[site].keys():
                continue
            site_cluster = Cluster(self.sites_data, self.intermediate, site)
            in_cluster = self.check_in_cluster(site_cluster)
            if in_cluster == True:
                sites_in_cluster.append(site)
        
        ground_energy = 10
        ground_site = ''
        for site in sites_in_cluster:
            if bool(self.sites_data[site]["converged"]) == True:
                if self.sites_data[site]["final_energy"] < ground_energy:
                    ground_site = site
                else:
                    continue
            elif bool(self.sites_data[site]["converged"]) == False:
                continue
        if ground_site == '':
            ground_site = None
        return ground_site
    
    def matching_site_clusters(self, site:Site, intermediate:str):
        site_clusters = site.unique_clusters(intermediate)
        matching_clusters = []
        for site_cluster in site_clusters:
            if self.check_in_cluster(site_cluster):
                matching_clusters.append(site_cluster)
        return matching_clusters

    def matching_ground_site(self, site:Site, intermediate:str):
        # this method will match the cluster object with a site object and return the ground state
        # for this cluster instance that matches the site passed into the method
        matching_clusters = self.matching_site_clusters(site, intermediate)
        print("matching clusters", matching_clusters)
        ground_sites = []
        for cluster in matching_clusters:
            ground_site = cluster.ground_site()
            print(cluster, ground_site)
            if cluster.ground_site() != None:
                ground_sites.append(ground_site)
            else:
                continue
        return ground_sites
