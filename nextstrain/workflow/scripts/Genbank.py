import time
import requests
import warnings
from dateutil import parser
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


class CDS:
    """
    A class to hold a coding sequence from a GenBank record.
    """
    def __init__(self, record: SeqRecord, feature: SeqFeature):
        """
        Initialize a CDS object from a GenBank record and feature.
        
        Parameters
        ----------
        record : SeqRecord
            The GenBank record containing the CDS feature.
        feature : SeqFeature
            The CDS feature to extract from the record.
            
        Raises
        ------
        ValueError
            If the feature is not a CDS, if it's not found in the record,
            or if required qualifiers are missing or have multiple values.
        """
        self.record = record
        self.protein_record = None
        self.feature = feature
        if feature.type != "CDS":
            raise ValueError("Feature is not a CDS")
        if feature not in record.features:
            raise ValueError("Feature not found in the record")   
        # Extract the nucleotide sequence from the record
        self.sequence = feature.extract(record.seq)
        # Extract the important information from the qualifiers
        self.qualifiers = feature.qualifiers
        for key in ['codon_start', 'product', 'protein_id', 'translation']:
            if key not in self.qualifiers:
                print(f"Warning: Qualifier {key} is missing from the CDS feature for {self.record.id}")
                value = ["Unknown"]
            else:
                value = self.qualifiers[key]
            if len(value) > 1:
                raise ValueError(f"Qualifier {key} has multiple values in the CDS feature for {self.record.id}")
            setattr(self, key, value[0])
    

    def __str__(self):
        """
        String representation of the CDS object.

        Returns:
        -------
        str
            A string representation of the CDS object.
        """
        return f"CDS for {self.product} ({self.protein_id}) for {self.record.id}\n{self.feature}"
        
    
    def __repr__(self):
        """
        Print default representation of the CDS object.
        """
        return f"CDS(record={self.record.id}, feature={self.protein_id}, product='{self.product}')"
    

    def __len__(self):
        """
        Length of the CDS nucleotide sequence.

        Returns:
        -------
        int
            Length of the nucleotide sequence.
        """
        return len(self.sequence)
    

    @property
    def ambiguous(self):
        """
        Count how many ambiguous nucleotides are in the sequence.

        Returns:
        -------
        int
            Number of ambiguous nucleotides in the sequence.
        """
        return sum(1 for base in self.sequence if base not in "ACGT")


    def fetch_protein_record(self, email="example@default.com"):
        """
        Fetch the GenBank record with protein id.

        Parameters
        ----------
        email : str
            Email address for NCBI Entrez API usage

        Returns:
        -------
        Bio.SeqRecord.SeqRecord or None
            The protein record if fetched successfully, otherwise None.
        """
        Entrez.email = email
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                handle = Entrez.efetch(db="protein", id=self.protein_id, rettype="gb", retmode="text")
                protein_record = SeqIO.read(handle, "genbank")
                for warning in w:
                    if "BiopythonParserWarning" in str(warning.category):
                        print(f"BiopythonParserWarning: Found malformed header in GenBank record for {self.protein_id}")
                handle.close()
            self.protein_record = protein_record
            return protein_record
        except Exception as e:
            print(f"Error fetching record for accession {self.protein_id}: {e}")
            return None
        
    
    def check_translation(self, print_differences=True):
        """
        Check if the translated sequence matches the translation qualifier.

        This isn't always the case, so this is a warning rather than an error.
        Sometimes there is read-through or non-standard translation.

        Parameters
        ----------
        print_differences : bool
            If True, print the differences between the sequences if they do not match.

        Returns:
        -------
        bool
            True if the translation matches, False otherwise.
        """
        # Check if there are translation-related qualifiers
        if 'transl_except' in self.qualifiers:
            transl_except = self.qualifiers.get('transl_except')
            print(f"Warning: There is a 'transl_except' qualifier ({transl_except[0]}), indicating non-standard translation.")
        # Check at the translation is the same as the translated sequence
        translated_seq = str(self.sequence[(int(self.codon_start) - 1):].translate(to_stop=False))
        if translated_seq[-1] == '*':
            translated_seq = translated_seq[:-1]
        if translated_seq != self.translation:
            print("Warning: Translated sequence does not match the translation qualifier:")
            # Check if the lengths are the same
            if len(translated_seq) != len(self.translation):
                print(f"Length of translated sequence ({len(translated_seq)}) does not match length of translation qualifier ({len(self.translation)}).")
            else:
                print("Sequences are the same length but differ in sequence.")
            # Check at which positions the sequences differ
            differences = []
            for i, (a, b) in enumerate(zip(translated_seq, self.translation)):
                if a != b:
                    differences.append((i + 1, f"{a} -> {b}"))
            if differences and print_differences:
                print("Differences at positions:")
                for pos, diff in differences:
                    print(f"  Position {pos}: {diff}")
            return False
        return True
    

class GenBankRecord:
    """
    A class for fetching and parsing GenBank records from NCBI.
    
    An interface to retrieve GenBank records using Biopython's
    Entrez module and offers convenient methods to access key
    metadata, annotations, and features within these records.
    """
    def __init__(self, accession, record=None, email="example@default.com"):
        """
        Initialize the GenBankRecord object with accession.

        Parameters
        ----------
        accession : str
            The GenBank accession number.
        record : Bio.SeqRecord.SeqRecord, optional
            A pre-fetched GenBank record object.
        email : str
            Email address for NCBI Entrez API usage
        
        Raises
        ------
        ValueError
            If the supplied record's accession does not match the provided accession.
        TypeError
            If the supplied record is not a Bio.SeqRecord.SeqRecord object.
        """
        self.accession = accession
        self.url = f"https://www.ncbi.nlm.nih.gov/nucleotide/{self.accession}"
        if record is not None:
            if record.id != accession:
                raise ValueError(f"Supplied record accession {record.id} does not match the provided accession {accession}.")
            if not isinstance(record, SeqRecord):
                raise TypeError(f"Expected Bio.SeqRecord.SeqRecord, got {type(record).__name__} instead")
        self.record = record
        self.email = email
        Entrez.email = self.email


    def __str__(self):
        """
        String representation of the GenBankRecord object.

        Returns:
        -------
        str
            A string representation of the GenBankRecord object.
        """
        if self.record is not None:
            return f"GenBankRecord(accession={self.record.id})"
        return f"GenBankRecord(accession={self.accession})"
    

    def __repr__(self):
        """
        Print default representation of the GenBankRecord object.
        """
        return f"GenBankRecord(accession={self.accession}, record={self.record})"


    def __len__(self):
        """
        Length of the nucleotide sequence in the GenBank record.

        Returns:
        -------
        int
            Length of the nucleotide sequence.
        """
        if self.record is not None:
            return len(self.record.seq)
        return 0
    

    def fetch(self, reties=3, retry_delay=5):
        """
        Fetch the GenBank record from NCBI using the accession number.

        Returns:
        -------
        GenBankRecord or None
            The GenBankRecord object if fetched successfully, otherwise None.
        """
        attempts = 0
        while attempts < reties:
            try:
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter("always")
                    handle = Entrez.efetch(db="nucleotide", id=self.accession, rettype="gb", retmode="text")
                    self.record = SeqIO.read(handle, "genbank")
                    for warning in w:
                        if "BiopythonParserWarning" in str(warning.category) and "malformed header line" in str(warning.message):
                            print(f"BiopythonParserWarning: Found malformed header in GenBank record {self.accession}")
                    handle.close()
                return self
            except Exception as e:
                attempts += 1
                if attempts < reties:
                    print(f"Error fetching record for accession {self.accession}: {e}. Retrying ({attempts}/{reties})...")
                    time.sleep(retry_delay)
                else:
                    print(f"Failed to fetch record for accession {self.accession} after {reties} attempts: {e}")
                    return None
            
    
    @property
    def source(self):
        """
        Parse the source from the record's features.

        Returns:
        -------
        dict or None
            A dictionary of source qualifiers if found, otherwise None.
        """
        if self.record is None:
            return None
        sources = []
        for feature in self.record.features:
            if feature.type == "source":
                sources.append(feature)
        if not sources:
            raise ValueError(f"No source feature found in record {self.record.id}.")
        elif len(sources) > 1:
            raise ValueError(f"Multiple source features found in record {self.record.id}.")
        else:
            source = sources[0]
            return source.qualifiers


    @property
    def references(self):
        """
        Parse the references from the record header.

        Returns:
        -------
        list or None
            A list of Bio.SeqFeature.Reference objects, otherwise None.
        """
        if self.record is None:
            return None
        if "references" in self.record.annotations:
            references = self.record.annotations["references"]
            return references
        else:
            raise ValueError(f"No references found in record {self.record.id}.")


    @property
    def authors(self):
        """
        Parse the authors from the record's reference information.
        
        There can be multiple references, this method chooses the 
        longest list of authors if multiple references are present.

        Returns:
        -------
        str or None
            A string of the authors if found, otherwise None.
        """
        if self.record is None:
            return None
        authors = None
        for reference in self.references:
            if hasattr(reference, 'authors'):
                if authors is None or len(reference.authors) > len(authors):
                    authors = reference.authors
        return authors
    

    @property
    def title(self):
        """
        Parse the title from the record's reference information.
        
        There can be multiple references, this method chooses looks
        for the common title "Direct Submission". If another title is found,
        it will return that title. Otherwise, it returns Direct Submission.

        Returns:
        -------
        str or None
            A string of the title if found, otherwise "Direct Submission".
        """
        if self.record is None:
            return None
        title = None
        for reference in self.references:
            if hasattr(reference, 'title'):
                if reference.title.lower() != "direct submission":
                    title = reference.title
                break
        if title is None:
            title = "Direct Submission"
        return title


    @property
    def journal(self):
        """
        Parse the journal from the record's reference information.

        Returns:
        -------
        str or None
            A string of the journal if found, otherwise Unpublished.
        """
        if self.record is None:
            return None
        journal = None
        for reference in self.references:
            if reference.title.lower() != "direct submission":
                if hasattr(reference, 'journal'):
                    journal = reference.journal
                break
        if journal is None:
            journal = "Unpublished"
        return journal


    @property
    def paper_link(self):
        """
        Parse the paper link from the record's reference information.

        Returns:
        -------
        str or None
            A string of the PubMed link if found, otherwise None.
        """
        if self.record is None:
            return None
        paper_link = None
        for reference in self.references:
            if reference.title.lower() != "direct submission":
                if hasattr(reference, 'pubmed_id') and reference.pubmed_id:
                    paper_link = f"https://pubmed.ncbi.nlm.nih.gov/{reference.pubmed_id}"
                break
        return paper_link


    @property
    def submission(self):
        """
        Parse the submission information from the record's reference.

        Returns:
        -------
        str or None
            A string of the journal if found, otherwise None.
        """
        if self.record is None:
            return None
        submission = None
        for reference in self.references:
            if reference.title.lower() == "direct submission":
                if hasattr(reference, 'journal'):
                    if reference.journal.startswith("Submitted"):
                        submission = reference.journal
                        break
        return submission
    

    @property
    def strain(self):
        """
        Parse the strain name from the GenBank record.

        This is usually found in the source feature.
        It's either the 'strain' or 'isolate' qualifier.
        This looks for both, but gives priority to 'strain'.

        Returns:
        -------
        str or None
            A string of the strain name if found, otherwise None.
        """
        if self.record is None:
            return None
        source = self.source
        if source is not None:
            if 'strain' in source:
                return source['strain'][0]
            elif 'isolate' in source:
                return source['isolate'][0]
        return None
    

    @property
    def organism(self):
        """
        Parse the virus name from the GenBank record.

        This is usually found in the source feature.
        But it can be found in the 'organism' annotation.
        Or in the 'source' annotation.

        Returns:
        -------
        str or None
            A string of the organism name if found, otherwise None.
        """
        if self.record is None:
            return None
        source = self.source
        if source is not None:
            if 'organism' in source:
                return source['organism'][0]
            elif 'organism' in self.record.annotations:
                return self.record.annotations['organism']
            elif 'source' in self.record.annotations:
                return self.record.annotations['source']
        return None
    

    @property
    def host(self):
        """
        Parse the host that the virus was isolated from.

        Returns:
        -------
        str or None
            A string of the host name if found, otherwise Unknown.
        """
        if self.record is None:
            return None
        source = self.source
        if 'host' in source:
            return source['host'][0]
        return "Unknown"
    

    @property
    def date(self):
        """
        Parse the collection date from the GenBank record.

        Returns:
        -------
        str or None
            A string of the data after some formatting
        """
        if self.record is None:
            return None
        source = self.source
        if 'collection_date' in source:
            unformatted_date =  source['collection_date'][0]
            try:
                formatted_date = parser.parse(unformatted_date)
                parts = max(
                    len(unformatted_date.split("-")),
                    len(unformatted_date.split("/")),
                    len(unformatted_date.split("_"))
                )
                if parts == 3:
                    return formatted_date.strftime("%Y-%m-%d")
                elif parts == 2:
                    return formatted_date.strftime("%Y-%m") + "-XX"
                else:
                    return formatted_date.strftime("%Y") + "-XX-XX"
            except Exception as e:
                print(f"Error parsing collection date for {self.accession}: {e}")
                return unformatted_date
        return None


    @property
    def location(self):
        """
        Parse the unformatted location from the GenBank record.

        Generally in the format `Country: Region`

        Returns:
        -------
        str or None
            A string of the location if found, otherwise None.
        """
        if self.record is None:
            return None
        source = self.source
        if "geo_loc_name" in source:
            return source['geo_loc_name'][0]
        return None
    

    @property
    def country(self):
        """
        Parse the country from the location.

        Returns:
        -------
        str or None
            A string of the country if found, otherwise Unknown.
        """
        if self.record is None:
            return None
        location = self.location
        if location is not None:
            parts = location.split(":")
            if len(parts) > 0:
                return parts[0].strip().lower().capitalize()
        return "Unknown"


    @property
    def local(self):
        """
        Parse the local area from the location.

        Returns:
        -------
        str
            A string of the local area if found, otherwise Unknown.
        """
        if self.record is None:
            return None
        location = self.location
        if location is not None:
            parts = location.split(":")
            if len(parts) > 1:
                return parts[-1].strip().lower().capitalize()
        return "Unknown"


    def fetch_geographic_information(self):
        """
        Fetch geographic information for the GenBank record.

        The country and local area are messy in their formatting.
        This method manages fuzzy matching by making a call to the
        [REST Countries API](https://restcountries.com/#rest-countries).
        

        Returns:
        -------
        dict
            A dictionary with geographic info from the REST API
        """
        if self.record is None:
            return None
        source = self.source
        if source is None or "geo_loc_name" not in source:
            return None
        
        location = source['geo_loc_name'][0].split(":")
        country = location[0].strip().lower().capitalize()
        if len(location) > 1:
            local = location[-1].strip().lower().capitalize()
        else:
            local = "Unknown"
        
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
            'Accept': 'application/json'
        }
        url = f"https://restcountries.com/v3.1/name/{country}"
        try:
            response = requests.get(url, headers=headers, timeout=10)
            if response.status_code == 200:
                data = response.json()
                if data:
                    return {
                        "country": data[0].get("name", {}).get("common", "Unknown"),
                        "region": data[0].get("region", "Unknown"),
                        "subregion": data[0].get("subregion", "Unknown"),
                        "local": local
                    }
        except Exception as e:
            print(f"Error fetching geographic information for {country}: {e}")
            return None
        
    
    @property
    def sequence(self):
        """
        Extract the nucleotide sequence from the GenBank record.

        Returns:
        -------
        Bio.Seq.Seq or None
            Seq object of the nucleotide sequence if found, otherwise None.
        """
        if self.record is None:
            return None
        return self.record.seq
    

    @property
    def ambiguous(self):
        """
        Count how many ambiguous nucleotides are in the sequence.
        """
        if self.record is None:
            return None
        return sum(1 for base in self.sequence if base not in "ACGT")
    

    @property
    def length(self):
        """
        Count how many nucleotides are in the sequence.
        """
        return self.__len__()


    @property
    def coding_regions(self):
        """
        Extract the coding regions from the GenBank record.
        """
        if self.record is None:
            return None
        coding_sequences = []
        for feature in self.record.features:
            if feature.type == "CDS":
                try:
                    coding_sequences.append(CDS(self.record, feature))
                except Exception as e:
                    print(f"Warning: Unable to parse CDS feature for {self.accession}:\n {e}")
                    continue
        return coding_sequences
