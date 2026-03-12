{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "187e09b7-5838-42c6-b6c2-bfa9efef2619",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Requirement already satisfied: pandas in c:\\programdata\\anaconda3\\lib\\site-packages (2.2.2)\n",
      "Requirement already satisfied: numpy in c:\\programdata\\anaconda3\\lib\\site-packages (1.26.4)\n",
      "Requirement already satisfied: matplotlib in c:\\programdata\\anaconda3\\lib\\site-packages (3.9.2)\n",
      "Requirement already satisfied: seaborn in c:\\programdata\\anaconda3\\lib\\site-packages (0.13.2)\n",
      "Requirement already satisfied: biopython in c:\\programdata\\anaconda3\\lib\\site-packages (1.86)\n",
      "Requirement already satisfied: requests in c:\\programdata\\anaconda3\\lib\\site-packages (2.32.3)\n",
      "Requirement already satisfied: python-dateutil>=2.8.2 in c:\\programdata\\anaconda3\\lib\\site-packages (from pandas) (2.9.0.post0)\n",
      "Requirement already satisfied: pytz>=2020.1 in c:\\programdata\\anaconda3\\lib\\site-packages (from pandas) (2024.1)\n",
      "Requirement already satisfied: tzdata>=2022.7 in c:\\programdata\\anaconda3\\lib\\site-packages (from pandas) (2023.3)\n",
      "Requirement already satisfied: contourpy>=1.0.1 in c:\\programdata\\anaconda3\\lib\\site-packages (from matplotlib) (1.2.0)\n",
      "Requirement already satisfied: cycler>=0.10 in c:\\programdata\\anaconda3\\lib\\site-packages (from matplotlib) (0.11.0)\n",
      "Requirement already satisfied: fonttools>=4.22.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from matplotlib) (4.51.0)\n",
      "Requirement already satisfied: kiwisolver>=1.3.1 in c:\\programdata\\anaconda3\\lib\\site-packages (from matplotlib) (1.4.4)\n",
      "Requirement already satisfied: packaging>=20.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from matplotlib) (24.1)\n",
      "Requirement already satisfied: pillow>=8 in c:\\programdata\\anaconda3\\lib\\site-packages (from matplotlib) (10.4.0)\n",
      "Requirement already satisfied: pyparsing>=2.3.1 in c:\\programdata\\anaconda3\\lib\\site-packages (from matplotlib) (3.1.2)\n",
      "Requirement already satisfied: charset-normalizer<4,>=2 in c:\\programdata\\anaconda3\\lib\\site-packages (from requests) (3.3.2)\n",
      "Requirement already satisfied: idna<4,>=2.5 in c:\\programdata\\anaconda3\\lib\\site-packages (from requests) (3.7)\n",
      "Requirement already satisfied: urllib3<3,>=1.21.1 in c:\\programdata\\anaconda3\\lib\\site-packages (from requests) (2.2.3)\n",
      "Requirement already satisfied: certifi>=2017.4.17 in c:\\programdata\\anaconda3\\lib\\site-packages (from requests) (2024.8.30)\n",
      "Requirement already satisfied: six>=1.5 in c:\\programdata\\anaconda3\\lib\\site-packages (from python-dateutil>=2.8.2->pandas) (1.16.0)\n"
     ]
    }
   ],
   "source": [
    "!pip install pandas numpy matplotlib seaborn biopython requests\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "9c9af618-3e55-47d5-9847-b7efbbf87fa8",
   "metadata": {},
   "outputs": [],
   "source": [
    "import pandas as pd\n",
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "import seaborn as sns\n",
    "from Bio import Entrez, SeqIO\n",
    "import requests"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "3c5c51a4-d241-42c0-b874-c551879ddc52",
   "metadata": {},
   "source": [
    "Huntington’s Disease is caused by mutation in HTT gene.\n",
    "Specifically due to CAG repeat expansion.\n",
    "\n",
    "Normal: 10–35 repeats  \n",
    "Disease: >36 repeats"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "8289a264-5215-4ae1-a35b-e8217322c45d",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "{'Count': '2', 'RetMax': '2', 'RetStart': '0', 'IdList': ['6532', '3064'], 'TranslationSet': [{'From': 'Homo sapiens[Organism]', 'To': '\"Homo sapiens\"[Organism]'}], 'TranslationStack': [{'Term': 'HTT[Gene]', 'Field': 'Gene', 'Count': '1121', 'Explode': 'N'}, {'Term': '\"Homo sapiens\"[Organism]', 'Field': 'Organism', 'Count': '360486', 'Explode': 'Y'}, 'AND'], 'QueryTranslation': 'HTT[Gene] AND \"Homo sapiens\"[Organism]'}"
      ]
     },
     "execution_count": 9,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "Entrez.email = \"priyas.162001@gmail.com\"\n",
    "\n",
    "handle = Entrez.esearch(db=\"gene\", term=\"HTT[Gene] AND Homo sapiens[Organism]\")\n",
    "record = Entrez.read(handle)\n",
    "\n",
    "record"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "id": "5281e085-5dd6-4ecc-bcaf-9681ad5cb2a3",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "'This gene encodes an integral membrane protein that transports the neurotransmitter serotonin from synaptic spaces into presynaptic neurons. The encoded protein terminates the action of serotonin and recycles it in a sodium-dependent manner. This protein is a target of psychomotor stimulants, such as amphetamines and cocaine, and is a member of the sodium:neurotransmitter symporter family. A repeat length polymorphism in the promoter of this gene has been shown to affect the rate of serotonin uptake. There have been conflicting results in the literature about the possible effect, if any, that this polymorphism may play in behavior and depression. [provided by RefSeq, May 2019]'"
      ]
     },
     "execution_count": 11,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "gene_id = record['IdList'][0]\n",
    "\n",
    "handle = Entrez.efetch(db=\"gene\", id=gene_id, retmode=\"xml\")\n",
    "gene_data = Entrez.read(handle)\n",
    "\n",
    "gene_data[0]['Entrezgene_summary']"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "id": "151b76d3-fc5f-45a6-ae41-7174677dd4e6",
   "metadata": {},
   "outputs": [],
   "source": [
    "gene_id = record['IdList'][1]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "id": "c797846f-daa4-43d3-bc4b-41d20ede0843",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "\"Huntingtin is a disease gene linked to Huntington's disease, a neurodegenerative disorder characterized by loss of striatal neurons. This is thought to be caused by an expanded, unstable trinucleotide repeat in the huntingtin gene, which translates as a polyglutamine repeat in the protein product. A fairly broad range of trinucleotide repeats (9-35) has been identified in normal controls, and repeat numbers in excess of 40 have been described as pathological. The huntingtin locus is large, spanning 180 kb and consisting of 67 exons. The huntingtin gene is widely expressed and is required for normal development. It is expressed as 2 alternatively polyadenylated forms displaying different relative abundance in various fetal and adult tissues. The larger transcript is approximately 13.7 kb and is expressed predominantly in adult and fetal brain whereas the smaller transcript of approximately 10.3 kb is more widely expressed. The genetic defect leading to Huntington's disease may not necessarily eliminate transcription, but may confer a new property on the mRNA or alter the function of the protein. One candidate is the huntingtin-associated protein-1, highly expressed in brain, which has increased affinity for huntingtin protein with expanded polyglutamine repeats. This gene contains an upstream open reading frame in the 5' UTR that inhibits expression of the huntingtin gene product through translational repression. [provided by RefSeq, Jul 2016]\""
      ]
     },
     "execution_count": 15,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "gene_id = record['IdList'][1]\n",
    "\n",
    "handle = Entrez.efetch(db=\"gene\", id=gene_id, retmode=\"xml\")\n",
    "gene_data = Entrez.read(handle)\n",
    "\n",
    "gene_data[0]['Entrezgene_summary']"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "id": "12fdfb4c-a276-4d66-a794-c519b76a36d2",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "NM_002111.8 Homo sapiens huntingtin (HTT), transcript variant 2, mRNA\n"
     ]
    }
   ],
   "source": [
    "from Bio import Entrez, SeqIO\n",
    "\n",
    "Entrez.email = \"your_email@gmail.com\"\n",
    "\n",
    "handle = Entrez.efetch(\n",
    "    db=\"nucleotide\",\n",
    "    id=\"NM_002111\",\n",
    "    rettype=\"fasta\",\n",
    "    retmode=\"text\"\n",
    ")\n",
    "\n",
    "seq_record = SeqIO.read(handle, \"fasta\")\n",
    "\n",
    "print(seq_record.description)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 19,
   "id": "dd4a73c1-a443-44c4-a138-999f8f10089c",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "13498"
      ]
     },
     "execution_count": 19,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "len(seq_record.seq)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 21,
   "id": "45041b83-6e8b-4e89-96e1-f6f0c7b6d25e",
   "metadata": {},
   "outputs": [],
   "source": [
    "sequence = str(seq_record.seq)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 23,
   "id": "3355a034-9402-439f-ba91-dba92666ca9d",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "28"
      ]
     },
     "execution_count": 23,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "import re\n",
    "\n",
    "cag_blocks = re.findall(r\"(?:CAG){2,}\", sequence)\n",
    "\n",
    "len(cag_blocks)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 25,
   "id": "88868557-4b15-4c2a-8c84-522797440a6d",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[21,\n",
       " 2,\n",
       " 3,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 3,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2,\n",
       " 2]"
      ]
     },
     "execution_count": 25,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "repeat_counts = [len(block)//3 for block in cag_blocks]\n",
    "\n",
    "repeat_counts"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 27,
   "id": "937378b7-4e7c-434c-bf5d-cbf100844a7d",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "21"
      ]
     },
     "execution_count": 27,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "max(repeat_counts)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 29,
   "id": "71385c7c-2cbe-41f4-8a60-1dd2951185d7",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[(196, 259, 21),\n",
       " (968, 974, 2),\n",
       " (1126, 1135, 3),\n",
       " (1291, 1297, 2),\n",
       " (1533, 1539, 2),\n",
       " (1616, 1622, 2),\n",
       " (2073, 2079, 2),\n",
       " (2490, 2496, 2),\n",
       " (2664, 2673, 3),\n",
       " (3170, 3176, 2),\n",
       " (3500, 3506, 2),\n",
       " (4927, 4933, 2),\n",
       " (5011, 5017, 2),\n",
       " (5537, 5543, 2),\n",
       " (5713, 5719, 2),\n",
       " (6691, 6697, 2),\n",
       " (7735, 7741, 2),\n",
       " (8481, 8487, 2),\n",
       " (8626, 8632, 2),\n",
       " (9163, 9169, 2),\n",
       " (9381, 9387, 2),\n",
       " (10396, 10402, 2),\n",
       " (10500, 10506, 2),\n",
       " (11385, 11391, 2),\n",
       " (11608, 11614, 2),\n",
       " (11622, 11628, 2),\n",
       " (12801, 12807, 2),\n",
       " (13329, 13335, 2)]"
      ]
     },
     "execution_count": 29,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "positions = []\n",
    "\n",
    "for match in re.finditer(r\"(?:CAG){2,}\", sequence):\n",
    "    start = match.start()\n",
    "    end = match.end()\n",
    "    length = (end - start)//3\n",
    "    positions.append((start, end, length))\n",
    "\n",
    "positions"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 31,
   "id": "8ebe0d80-8539-4f67-a268-2431d4818755",
   "metadata": {},
   "outputs": [],
   "source": [
    "mutant_sequence = (\n",
    "    sequence[:196] +\n",
    "    \"CAG\"*45 +\n",
    "    sequence[259:]\n",
    ")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 33,
   "id": "90fe12aa-c445-4498-bbd0-5f809f9473e3",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "45"
      ]
     },
     "execution_count": 33,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "mutant_blocks = re.findall(r\"(?:CAG){2,}\", mutant_sequence)\n",
    "\n",
    "mutant_counts = [len(b)//3 for b in mutant_blocks]\n",
    "\n",
    "max(mutant_counts)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 35,
   "id": "1165da87-a49b-454a-895e-b95bdeb6303c",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>Sequence Type</th>\n",
       "      <th>Max CAG Repeat</th>\n",
       "      <th>Disease Status</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>Normal HTT</td>\n",
       "      <td>21</td>\n",
       "      <td>Normal (&lt;27)</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>Mutant HTT</td>\n",
       "      <td>45</td>\n",
       "      <td>Huntington Risk (&gt;36)</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "  Sequence Type  Max CAG Repeat         Disease Status\n",
       "0    Normal HTT              21           Normal (<27)\n",
       "1    Mutant HTT              45  Huntington Risk (>36)"
      ]
     },
     "execution_count": 35,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "normal_max = max(repeat_counts)\n",
    "mutant_max = max(mutant_counts)\n",
    "\n",
    "comparison = {\n",
    "    \"Sequence Type\": [\"Normal HTT\", \"Mutant HTT\"],\n",
    "    \"Max CAG Repeat\": [normal_max, mutant_max],\n",
    "    \"Disease Status\": [\n",
    "        \"Normal (<27)\",\n",
    "        \"Huntington Risk (>36)\"\n",
    "    ]\n",
    "}\n",
    "\n",
    "import pandas as pd\n",
    "\n",
    "df_compare = pd.DataFrame(comparison)\n",
    "\n",
    "df_compare"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 37,
   "id": "46aa01a7-92ae-4170-96d9-2800eec7f9c9",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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nj179PzzzysiIkKlSpWyr5RHjRqlVq1aacGCBdq/f7+effZZ3X///fLw8FD9+vX1/vvva9euXRo1apQCAgI0ffp0+7z79u1T79697UHw22+/1csvv6yffvpJ8+bNc+m19+vXT6+//rqWLl2qRx55xD5+8uRJffrppxo4cKA8PT21f/9+derUSXfddZfmzZunMmXK6I8//tDq1auVlZUlPz+/yz5Xr169NGTIEL3//vvq37+/pAu7s3r06HFVuxq+/vprtW7dWq1atbL/ssqd7/jx45KkMWPGKCwsTKdPn9Ynn3yi6OhorV271h7Wcs2aNUs333yzpk2bJkl64YUX1LFjRyUlJal06dJ64YUXlJ6ero8++khff/21/XHly5e/bJ/nz5+/4tdY1Hbs2KHjx4/r0UcfdQi0lzN37lxJFz434eHhGjx4sObOnesQqL766itJ0j333FO0TV8kKSlJZcqUUXBwsH3s8ccf14IFCzRo0CBNmjRJx48f17hx49SsWTN9++23Cg0NtdcmJyerV69eeu655zRu3Dj77twTJ05o5syZki78ou/atas2b96s4cOHq1mzZvrtt980ZswYRUdHa8eOHfL19ZUk/fzzz+rYsaP9PzU//fSTJk2apO3bt2vdunWS5PJnx9WfUUZGhlq1aqV9+/Zp7NixqlevnjZv3qz4+HglJiY6hY+VK1dq165dmjhxomw2m0aMGKFOnTopNjZWv/76q2bOnKnU1FQNHTpU9913nxITEx36KMx7KF34D87jjz+uypUrS5K2bdump556Sn/88YdefPFFh57yWmfmZfLkyYqLi9Pzzz+vFi1a6Ny5c/rpp58cQvnGjRsVExOjevXqae7cufL29tbs2bPVpUsXvf/+++rZs6fDnI888og6deqk9957TwcPHtSwYcP04IMP2n9+KCLWbojD9eTiXTiZmZkmMjLSREVFmZycHGOM866Q1atXG0lm8uTJDvN88MEHRpJ566237GNVqlQx7u7uZu/evQ6169evN5JMly5dHMYHDx5sJJlBgwY5jN97773mpptuyvc1ZGdnm3Pnzpl33nnHuLu7m+PHj9uXFXZ3ZKNGjUyzZs0cxmbPnm0kme+++84YY8xHH31kJDntsi2Mi9/H2NhYExUVZYwx5ocffjCSzIYNG0xCQoLT7sj8+h8zZozTLp3C7o48f/68OXfunGnTpo3p1q2bfTwpKclIMrfeeqs5f/68fXz79u1Gknn//fftY1eyO1JSgbdrvTtyyZIlRpJ54403Cv060tPTTWBgoLn99tsdXpvNZjO//PKLfeyJJ54wksxPP/3k8PicnBxz7tw5++3i9zk/uZ+d3MccPnzYvPjii069f/3110aSefXVVx0ef/DgQePr62uGDx/uMKck8+mnnzrUPvroo8bNzc389ttvxhhj3n//fSPJfPzxxw51uZ/V2bNn59lz7uvcuHGjkWS+/fZb+zJXPjuu/ozeeOMNI8ksXbrUYTx3l+aXX35pH5NkwsLCzOnTp+1jy5cvN5JMgwYN7OtAY4yZNm2akWR2795tHyvse3ip3PXVuHHjTFBQkMPz5LfOzF128Xekc+fOpkGDBgW+H7fffrsJCQkxp06dso+dP3/e1K1b11SqVMn+3LnfnQEDBjg8fvLkyUaSOXz4cIHPA9ewOxJ58vLy0vjx47Vjxw4tXbo0z5rc/xFdepbO//zP/6hUqVJau3atw3i9evVUs2bNPOfq3Lmzw/3atWtLkjp16uQ0fvz4cYddkrt27dI999yjoKAgubu7y9PTUw8//LCys7P13//+9/Iv9hJ9+/bV1q1btXfvXvvY/Pnz1bhxY/vZoQ0aNJCXl5cee+wxLVy40GkXT2H169dPO3bs0Hfffae5c+eqWrVq+e6yKipvvPGGGjVqJB8fH3l4eMjT01Nr167Vnj17nGo7deokd3d3+/169epJUp5nz7rC19dXCQkJed5yt6Zc75YuXaq0tDT169fPPtavXz8ZYzR//vzLPv7TTz+Vp6en/Va6dOlCPe8PP/xgf0z58uU1btw4jRw5Uo8//ri9ZuXKlbLZbHrwwQd1/vx5+y0sLEz169fXhg0bHOYMCAhw2lLXu3dv5eTkaNOmTfY5y5Qpoy5dujjM2aBBA4WFhTnM+euvv6p3794KCwuzfydbtmwpSXl+zorDunXrVKpUKfXo0cNhPHd9den6qVWrVg5bmnLXQXfffbfDFq/c8Uu/A4V5D3P7atu2rUqXLm1/b1588UUdO3ZMKSkpDo8vaJ15sSZNmujbb7/VgAED9MUXXzgdTpGenq7//Oc/6tGjh8OZye7u7nrooYf0+++/O6zvJOctt0X13YcjQhjy1atXLzVq1EijR4/O8/iqY8eOycPDw2EXiHTh+IqwsDAdO3bMYbygXVQ33XSTw/3cM93yG8/IyJAkHThwQHfddZf++OMPvf7669q8ebMSEhLsx5qcPXu2MC/VwQMPPCBvb2/75SF+/PFHJSQkqG/fvvaaatWqac2aNQoJCdHAgQNVrVo1VatWzekYjMtp0aKFatSooTfffFOLFi1Sv379XNod5qqpU6eqf//+atq0qT7++GNt27ZNCQkJ6tChQ57vVVBQkMP93N3LV/K+XszNzU1RUVF53tzcrv1qKXfXUFJSUqEfM3fuXPn4+KhDhw46efKkTp48qXr16qlq1apasGCB/Vig3Lkv/eUVHR1tD56X/iekINWqVVNCQoK2b9+uDz/8UPXr11d8fLz9eEJJOnLkiIwxCg0NdQh6np6e2rZtm44ePeow58W7JnOFhYVJkv17fOTIEZ08eVJeXl5OcyYnJ9vnPH36tO666y795z//0fjx47VhwwYlJCRo2bJlkq78s+Pqz+jYsWMKCwtz+j6FhITIw8PDaf10peugXIV5D7dv36527dpJkt5++239+9//VkJCgkaPHi3J+b0pzG59SRo5cqSmTJmibdu26e6771ZQUJDatGljv2TKiRMnZIzJc74KFSo49JiruL77cMQxYciXzWbTpEmTFBMTk+ep70FBQTp//rz+/PNPhyBmjFFycrIaN27sNF9RW758udLT07Vs2TJVqVLFPp6YmHjFc5YtW1Zdu3bVO++8Y7/MgY+Pj+6//36Hurvuukt33XWXsrOztWPHDs2YMUODBw9WaGhongfW56dv3756/vnnZbPZFBsbm2+dj4+P0wHFkpx+oRZk8eLFio6O1pw5cxzGC3ONqpIsKipKN910kz799FPFx8df9rP63//+V1u2bJH0Vzi41BdffKGOHTsqJiZGo0aN0ooVK+y/gCWpTJky9stZXPoLryA+Pj72xzVu3FitWrXSLbfcosGDB6tz587y9/dXuXLlZLPZtHnzZofjMnNdOnbkyBGnmuTkZIfecg/QXr16dZ59BQQESLqwpefQoUPasGGDfeuXJKeTdVzl6s8oKChI//nPf2SMcahNSUnR+fPnVa5cuavq51KFeQ+XLFkiT09PrVy5Uj4+Pva6/K7pV9h1poeHh4YOHaqhQ4fq5MmTWrNmjUaNGqX27dvr4MGDKlu2rNzc3HT48GGnx+YebF/U7wcKhy1hKFDbtm0VExOjcePGOZ2VmHt9q8WLFzuMf/zxx0pPT3e6/lVxyF1JXfxLxRijt99++6rm7du3rw4dOqRVq1Zp8eLF6tatm8qUKZNnrbu7u5o2bWrf+vbNN9+49FyxsbHq0qWLhg0bpooVK+ZbV7VqVaWkpDis7LOysvTFF1841Xp7e+f5P1abzeb0C3j37t0OB0a7qiT8D9nT01MjRozQTz/9pJdeeinPmpSUFP373/+W9NcB+W+//bbWr1/vcFu1apU8PT3tJ4VERUWpXbt2evvtt53OQi0KQUFBmjhxoo4cOaIZM2ZIurB73xijP/74I8+tjbfeeqvDHKdOndKKFSscxt577z25ubnZd4937txZx44dU3Z2dp5z1qpVS1Le30lJevPNN516d+Wz4+rPqE2bNjp9+rRTwHnnnXfsy4tSYd5Dm80mDw8Ph138Z8+eLdJrMpYpU0Y9evTQwIEDdfz4ce3fv1+lSpVS06ZNtWzZMof3OicnR4sXL1alSpUKtdsTRY8tYbisSZMm6bbbblNKSopuueUW+3hMTIzat2+vESNGKC0tTc2bN7efHdmwYUM99NBDxd5bTEyMvLy8dP/992v48OHKyMjQnDlzdOLEiauat127dqpUqZIGDBig5ORkh12R0oXjqtatW6dOnTqpcuXKysjIsP/Sbdu2rUvPVaFChUJd3b5nz5568cUX1atXLw0bNkwZGRmaPn260ynwknTrrbdqw4YN+te//qXy5csrICBAtWrVUufOnfXSSy9pzJgxatmypfbu3atx48YpIiLiis9WzP2FPmnSJN19991yd3dXvXr1HC6e+ncwbNgw7dmzR2PGjNH27dvVu3dv+4VAN23apLfeest+eYF33nlHtWvXdjiD9mJdunTRihUr7FuJFy9erPbt26tt27bq06eP2rdvr5CQEKWlpWn37t1as2bNVZ0R+/DDD2vq1KmaMmWKBg4cqObNm+uxxx5T3759tWPHDrVo0UKlSpXS4cOHtWXLFt166632M3KlC0Guf//+OnDggGrWrKlVq1bp7bffVv/+/e1b+nr16qV3331XHTt21NNPP60mTZrI09NTv//+u9avX6+uXbuqW7duatasmcqWLasnnnhCY8aMkaenp9599119++23Tn27+tkp7M+oefPmevjhhzVr1izFxsZq//79uvXWW7VlyxZNmDBBHTt2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",
      "text/plain": [
       "<Figure size 700x500 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "import matplotlib.pyplot as plt\n",
    "\n",
    "plt.figure(figsize=(7,5))\n",
    "plt.bar(df_compare[\"Sequence Type\"], df_compare[\"Max CAG Repeat\"])\n",
    "plt.axhline(y=36, color='red', linestyle='--', label='Disease Threshold')\n",
    "plt.title(\"Normal vs Mutant HTT CAG Repeat Comparison\")\n",
    "plt.ylabel(\"Repeat Count\")\n",
    "plt.legend()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "a786c33a-1506-4067-a98e-9253291b0001",
   "metadata": {},
   "source": [
    "Comparison shows mutant HTT allele crosses pathogenic CAG threshold.\n",
    "This demonstrates computational modeling of Huntington disease mutation."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 42,
   "id": "69dec534-ea6d-4365-9a2f-d8f4e90bb7cb",
   "metadata": {},
   "outputs": [],
   "source": [
    "def huntington_risk_predictor(cag_repeat):\n",
    "    \n",
    "    if cag_repeat < 27:\n",
    "        return \"Normal\"\n",
    "    \n",
    "    elif cag_repeat < 36:\n",
    "        return \"Intermediate Risk\"\n",
    "    \n",
    "    else:\n",
    "        return \"Huntington Disease Risk\""
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 44,
   "id": "f38b4ccf-eacd-46af-84ad-11be0c4c45a7",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "'Normal'"
      ]
     },
     "execution_count": 44,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "huntington_risk_predictor(21)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 46,
   "id": "7dd5c839-12c4-4818-87f0-8b1c756b053c",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "'Huntington Disease Risk'"
      ]
     },
     "execution_count": 46,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "huntington_risk_predictor(45)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 48,
   "id": "32a13a5a-d29e-4020-85d7-aba19c042e34",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>Sequence Type</th>\n",
       "      <th>Max CAG Repeat</th>\n",
       "      <th>Disease Status</th>\n",
       "      <th>Predicted Status</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>Normal HTT</td>\n",
       "      <td>21</td>\n",
       "      <td>Normal (&lt;27)</td>\n",
       "      <td>Normal</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>Mutant HTT</td>\n",
       "      <td>45</td>\n",
       "      <td>Huntington Risk (&gt;36)</td>\n",
       "      <td>Huntington Disease Risk</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "  Sequence Type  Max CAG Repeat         Disease Status  \\\n",
       "0    Normal HTT              21           Normal (<27)   \n",
       "1    Mutant HTT              45  Huntington Risk (>36)   \n",
       "\n",
       "          Predicted Status  \n",
       "0                   Normal  \n",
       "1  Huntington Disease Risk  "
      ]
     },
     "execution_count": 48,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "df_compare[\"Predicted Status\"] = df_compare[\"Max CAG Repeat\"].apply(huntington_risk_predictor)\n",
    "\n",
    "df_compare"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 50,
   "id": "ab6c3ecc-b483-4160-b6a3-2635ce12c745",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "1"
      ]
     },
     "execution_count": 50,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "from Bio.Seq import Seq\n",
    "from Bio.SeqRecord import SeqRecord\n",
    "from Bio import SeqIO\n",
    "\n",
    "mut_record = SeqRecord(\n",
    "    Seq(mutant_sequence),\n",
    "    id=\"HTT_MUTANT\",\n",
    "    description=\"Simulated Huntington Disease Allele\"\n",
    ")\n",
    "\n",
    "SeqIO.write(mut_record, \"HTT_mutant.fasta\", \"fasta\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 52,
   "id": "b612dafb-6a27-4e2b-80c9-283f45a73b10",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>Patient_ID</th>\n",
       "      <th>CAG_Repeats</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>1</td>\n",
       "      <td>53</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>2</td>\n",
       "      <td>43</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>2</th>\n",
       "      <td>3</td>\n",
       "      <td>29</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>3</th>\n",
       "      <td>4</td>\n",
       "      <td>22</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>4</th>\n",
       "      <td>5</td>\n",
       "      <td>35</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "   Patient_ID  CAG_Repeats\n",
       "0           1           53\n",
       "1           2           43\n",
       "2           3           29\n",
       "3           4           22\n",
       "4           5           35"
      ]
     },
     "execution_count": 52,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "import numpy as np\n",
    "import pandas as pd\n",
    "\n",
    "np.random.seed(42)\n",
    "\n",
    "patients = pd.DataFrame({\n",
    "    \"Patient_ID\": range(1,201),\n",
    "    \"CAG_Repeats\": np.random.randint(15, 55, 200)\n",
    "})\n",
    "\n",
    "patients.head()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 54,
   "id": "2301c869-0a2b-4703-8ac9-48dda1d422d9",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>Patient_ID</th>\n",
       "      <th>CAG_Repeats</th>\n",
       "      <th>Disease</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>1</td>\n",
       "      <td>53</td>\n",
       "      <td>1</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>2</td>\n",
       "      <td>43</td>\n",
       "      <td>1</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>2</th>\n",
       "      <td>3</td>\n",
       "      <td>29</td>\n",
       "      <td>0</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>3</th>\n",
       "      <td>4</td>\n",
       "      <td>22</td>\n",
       "      <td>0</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>4</th>\n",
       "      <td>5</td>\n",
       "      <td>35</td>\n",
       "      <td>0</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "   Patient_ID  CAG_Repeats  Disease\n",
       "0           1           53        1\n",
       "1           2           43        1\n",
       "2           3           29        0\n",
       "3           4           22        0\n",
       "4           5           35        0"
      ]
     },
     "execution_count": 54,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "def label_disease(x):\n",
    "    if x >= 36:\n",
    "        return 1\n",
    "    else:\n",
    "        return 0\n",
    "\n",
    "patients[\"Disease\"] = patients[\"CAG_Repeats\"].apply(label_disease)\n",
    "\n",
    "patients.head()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 56,
   "id": "61eb7ca3-8c2a-4893-b54a-91068d24d14d",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<style>#sk-container-id-1 {\n",
       "  /* Definition of color scheme common for light and dark mode */\n",
       "  --sklearn-color-text: black;\n",
       "  --sklearn-color-line: gray;\n",
       "  /* Definition of color scheme for unfitted estimators */\n",
       "  --sklearn-color-unfitted-level-0: #fff5e6;\n",
       "  --sklearn-color-unfitted-level-1: #f6e4d2;\n",
       "  --sklearn-color-unfitted-level-2: #ffe0b3;\n",
       "  --sklearn-color-unfitted-level-3: chocolate;\n",
       "  /* Definition of color scheme for fitted estimators */\n",
       "  --sklearn-color-fitted-level-0: #f0f8ff;\n",
       "  --sklearn-color-fitted-level-1: #d4ebff;\n",
       "  --sklearn-color-fitted-level-2: #b3dbfd;\n",
       "  --sklearn-color-fitted-level-3: cornflowerblue;\n",
       "\n",
       "  /* Specific color for light theme */\n",
       "  --sklearn-color-text-on-default-background: var(--sg-text-color, var(--theme-code-foreground, var(--jp-content-font-color1, black)));\n",
       "  --sklearn-color-background: var(--sg-background-color, var(--theme-background, var(--jp-layout-color0, white)));\n",
       "  --sklearn-color-border-box: var(--sg-text-color, var(--theme-code-foreground, var(--jp-content-font-color1, black)));\n",
       "  --sklearn-color-icon: #696969;\n",
       "\n",
       "  @media (prefers-color-scheme: dark) {\n",
       "    /* Redefinition of color scheme for dark theme */\n",
       "    --sklearn-color-text-on-default-background: var(--sg-text-color, var(--theme-code-foreground, var(--jp-content-font-color1, white)));\n",
       "    --sklearn-color-background: var(--sg-background-color, var(--theme-background, var(--jp-layout-color0, #111)));\n",
       "    --sklearn-color-border-box: var(--sg-text-color, var(--theme-code-foreground, var(--jp-content-font-color1, white)));\n",
       "    --sklearn-color-icon: #878787;\n",
       "  }\n",
       "}\n",
       "\n",
       "#sk-container-id-1 {\n",
       "  color: var(--sklearn-color-text);\n",
       "}\n",
       "\n",
       "#sk-container-id-1 pre {\n",
       "  padding: 0;\n",
       "}\n",
       "\n",
       "#sk-container-id-1 input.sk-hidden--visually {\n",
       "  border: 0;\n",
       "  clip: rect(1px 1px 1px 1px);\n",
       "  clip: rect(1px, 1px, 1px, 1px);\n",
       "  height: 1px;\n",
       "  margin: -1px;\n",
       "  overflow: hidden;\n",
       "  padding: 0;\n",
       "  position: absolute;\n",
       "  width: 1px;\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-dashed-wrapped {\n",
       "  border: 1px dashed var(--sklearn-color-line);\n",
       "  margin: 0 0.4em 0.5em 0.4em;\n",
       "  box-sizing: border-box;\n",
       "  padding-bottom: 0.4em;\n",
       "  background-color: var(--sklearn-color-background);\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-container {\n",
       "  /* jupyter's `normalize.less` sets `[hidden] { display: none; }`\n",
       "     but bootstrap.min.css set `[hidden] { display: none !important; }`\n",
       "     so we also need the `!important` here to be able to override the\n",
       "     default hidden behavior on the sphinx rendered scikit-learn.org.\n",
       "     See: https://github.com/scikit-learn/scikit-learn/issues/21755 */\n",
       "  display: inline-block !important;\n",
       "  position: relative;\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-text-repr-fallback {\n",
       "  display: none;\n",
       "}\n",
       "\n",
       "div.sk-parallel-item,\n",
       "div.sk-serial,\n",
       "div.sk-item {\n",
       "  /* draw centered vertical line to link estimators */\n",
       "  background-image: linear-gradient(var(--sklearn-color-text-on-default-background), var(--sklearn-color-text-on-default-background));\n",
       "  background-size: 2px 100%;\n",
       "  background-repeat: no-repeat;\n",
       "  background-position: center center;\n",
       "}\n",
       "\n",
       "/* Parallel-specific style estimator block */\n",
       "\n",
       "#sk-container-id-1 div.sk-parallel-item::after {\n",
       "  content: \"\";\n",
       "  width: 100%;\n",
       "  border-bottom: 2px solid var(--sklearn-color-text-on-default-background);\n",
       "  flex-grow: 1;\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-parallel {\n",
       "  display: flex;\n",
       "  align-items: stretch;\n",
       "  justify-content: center;\n",
       "  background-color: var(--sklearn-color-background);\n",
       "  position: relative;\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-parallel-item {\n",
       "  display: flex;\n",
       "  flex-direction: column;\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-parallel-item:first-child::after {\n",
       "  align-self: flex-end;\n",
       "  width: 50%;\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-parallel-item:last-child::after {\n",
       "  align-self: flex-start;\n",
       "  width: 50%;\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-parallel-item:only-child::after {\n",
       "  width: 0;\n",
       "}\n",
       "\n",
       "/* Serial-specific style estimator block */\n",
       "\n",
       "#sk-container-id-1 div.sk-serial {\n",
       "  display: flex;\n",
       "  flex-direction: column;\n",
       "  align-items: center;\n",
       "  background-color: var(--sklearn-color-background);\n",
       "  padding-right: 1em;\n",
       "  padding-left: 1em;\n",
       "}\n",
       "\n",
       "\n",
       "/* Toggleable style: style used for estimator/Pipeline/ColumnTransformer box that is\n",
       "clickable and can be expanded/collapsed.\n",
       "- Pipeline and ColumnTransformer use this feature and define the default style\n",
       "- Estimators will overwrite some part of the style using the `sk-estimator` class\n",
       "*/\n",
       "\n",
       "/* Pipeline and ColumnTransformer style (default) */\n",
       "\n",
       "#sk-container-id-1 div.sk-toggleable {\n",
       "  /* Default theme specific background. It is overwritten whether we have a\n",
       "  specific estimator or a Pipeline/ColumnTransformer */\n",
       "  background-color: var(--sklearn-color-background);\n",
       "}\n",
       "\n",
       "/* Toggleable label */\n",
       "#sk-container-id-1 label.sk-toggleable__label {\n",
       "  cursor: pointer;\n",
       "  display: block;\n",
       "  width: 100%;\n",
       "  margin-bottom: 0;\n",
       "  padding: 0.5em;\n",
       "  box-sizing: border-box;\n",
       "  text-align: center;\n",
       "}\n",
       "\n",
       "#sk-container-id-1 label.sk-toggleable__label-arrow:before {\n",
       "  /* Arrow on the left of the label */\n",
       "  content: \"▸\";\n",
       "  float: left;\n",
       "  margin-right: 0.25em;\n",
       "  color: var(--sklearn-color-icon);\n",
       "}\n",
       "\n",
       "#sk-container-id-1 label.sk-toggleable__label-arrow:hover:before {\n",
       "  color: var(--sklearn-color-text);\n",
       "}\n",
       "\n",
       "/* Toggleable content - dropdown */\n",
       "\n",
       "#sk-container-id-1 div.sk-toggleable__content {\n",
       "  max-height: 0;\n",
       "  max-width: 0;\n",
       "  overflow: hidden;\n",
       "  text-align: left;\n",
       "  /* unfitted */\n",
       "  background-color: var(--sklearn-color-unfitted-level-0);\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-toggleable__content.fitted {\n",
       "  /* fitted */\n",
       "  background-color: var(--sklearn-color-fitted-level-0);\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-toggleable__content pre {\n",
       "  margin: 0.2em;\n",
       "  border-radius: 0.25em;\n",
       "  color: var(--sklearn-color-text);\n",
       "  /* unfitted */\n",
       "  background-color: var(--sklearn-color-unfitted-level-0);\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-toggleable__content.fitted pre {\n",
       "  /* unfitted */\n",
       "  background-color: var(--sklearn-color-fitted-level-0);\n",
       "}\n",
       "\n",
       "#sk-container-id-1 input.sk-toggleable__control:checked~div.sk-toggleable__content {\n",
       "  /* Expand drop-down */\n",
       "  max-height: 200px;\n",
       "  max-width: 100%;\n",
       "  overflow: auto;\n",
       "}\n",
       "\n",
       "#sk-container-id-1 input.sk-toggleable__control:checked~label.sk-toggleable__label-arrow:before {\n",
       "  content: \"▾\";\n",
       "}\n",
       "\n",
       "/* Pipeline/ColumnTransformer-specific style */\n",
       "\n",
       "#sk-container-id-1 div.sk-label input.sk-toggleable__control:checked~label.sk-toggleable__label {\n",
       "  color: var(--sklearn-color-text);\n",
       "  background-color: var(--sklearn-color-unfitted-level-2);\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-label.fitted input.sk-toggleable__control:checked~label.sk-toggleable__label {\n",
       "  background-color: var(--sklearn-color-fitted-level-2);\n",
       "}\n",
       "\n",
       "/* Estimator-specific style */\n",
       "\n",
       "/* Colorize estimator box */\n",
       "#sk-container-id-1 div.sk-estimator input.sk-toggleable__control:checked~label.sk-toggleable__label {\n",
       "  /* unfitted */\n",
       "  background-color: var(--sklearn-color-unfitted-level-2);\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-estimator.fitted input.sk-toggleable__control:checked~label.sk-toggleable__label {\n",
       "  /* fitted */\n",
       "  background-color: var(--sklearn-color-fitted-level-2);\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-label label.sk-toggleable__label,\n",
       "#sk-container-id-1 div.sk-label label {\n",
       "  /* The background is the default theme color */\n",
       "  color: var(--sklearn-color-text-on-default-background);\n",
       "}\n",
       "\n",
       "/* On hover, darken the color of the background */\n",
       "#sk-container-id-1 div.sk-label:hover label.sk-toggleable__label {\n",
       "  color: var(--sklearn-color-text);\n",
       "  background-color: var(--sklearn-color-unfitted-level-2);\n",
       "}\n",
       "\n",
       "/* Label box, darken color on hover, fitted */\n",
       "#sk-container-id-1 div.sk-label.fitted:hover label.sk-toggleable__label.fitted {\n",
       "  color: var(--sklearn-color-text);\n",
       "  background-color: var(--sklearn-color-fitted-level-2);\n",
       "}\n",
       "\n",
       "/* Estimator label */\n",
       "\n",
       "#sk-container-id-1 div.sk-label label {\n",
       "  font-family: monospace;\n",
       "  font-weight: bold;\n",
       "  display: inline-block;\n",
       "  line-height: 1.2em;\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-label-container {\n",
       "  text-align: center;\n",
       "}\n",
       "\n",
       "/* Estimator-specific */\n",
       "#sk-container-id-1 div.sk-estimator {\n",
       "  font-family: monospace;\n",
       "  border: 1px dotted var(--sklearn-color-border-box);\n",
       "  border-radius: 0.25em;\n",
       "  box-sizing: border-box;\n",
       "  margin-bottom: 0.5em;\n",
       "  /* unfitted */\n",
       "  background-color: var(--sklearn-color-unfitted-level-0);\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-estimator.fitted {\n",
       "  /* fitted */\n",
       "  background-color: var(--sklearn-color-fitted-level-0);\n",
       "}\n",
       "\n",
       "/* on hover */\n",
       "#sk-container-id-1 div.sk-estimator:hover {\n",
       "  /* unfitted */\n",
       "  background-color: var(--sklearn-color-unfitted-level-2);\n",
       "}\n",
       "\n",
       "#sk-container-id-1 div.sk-estimator.fitted:hover {\n",
       "  /* fitted */\n",
       "  background-color: var(--sklearn-color-fitted-level-2);\n",
       "}\n",
       "\n",
       "/* Specification for estimator info (e.g. \"i\" and \"?\") */\n",
       "\n",
       "/* Common style for \"i\" and \"?\" */\n",
       "\n",
       ".sk-estimator-doc-link,\n",
       "a:link.sk-estimator-doc-link,\n",
       "a:visited.sk-estimator-doc-link {\n",
       "  float: right;\n",
       "  font-size: smaller;\n",
       "  line-height: 1em;\n",
       "  font-family: monospace;\n",
       "  background-color: var(--sklearn-color-background);\n",
       "  border-radius: 1em;\n",
       "  height: 1em;\n",
       "  width: 1em;\n",
       "  text-decoration: none !important;\n",
       "  margin-left: 1ex;\n",
       "  /* unfitted */\n",
       "  border: var(--sklearn-color-unfitted-level-1) 1pt solid;\n",
       "  color: var(--sklearn-color-unfitted-level-1);\n",
       "}\n",
       "\n",
       ".sk-estimator-doc-link.fitted,\n",
       "a:link.sk-estimator-doc-link.fitted,\n",
       "a:visited.sk-estimator-doc-link.fitted {\n",
       "  /* fitted */\n",
       "  border: var(--sklearn-color-fitted-level-1) 1pt solid;\n",
       "  color: var(--sklearn-color-fitted-level-1);\n",
       "}\n",
       "\n",
       "/* On hover */\n",
       "div.sk-estimator:hover .sk-estimator-doc-link:hover,\n",
       ".sk-estimator-doc-link:hover,\n",
       "div.sk-label-container:hover .sk-estimator-doc-link:hover,\n",
       ".sk-estimator-doc-link:hover {\n",
       "  /* unfitted */\n",
       "  background-color: var(--sklearn-color-unfitted-level-3);\n",
       "  color: var(--sklearn-color-background);\n",
       "  text-decoration: none;\n",
       "}\n",
       "\n",
       "div.sk-estimator.fitted:hover .sk-estimator-doc-link.fitted:hover,\n",
       ".sk-estimator-doc-link.fitted:hover,\n",
       "div.sk-label-container:hover .sk-estimator-doc-link.fitted:hover,\n",
       ".sk-estimator-doc-link.fitted:hover {\n",
       "  /* fitted */\n",
       "  background-color: var(--sklearn-color-fitted-level-3);\n",
       "  color: var(--sklearn-color-background);\n",
       "  text-decoration: none;\n",
       "}\n",
       "\n",
       "/* Span, style for the box shown on hovering the info icon */\n",
       ".sk-estimator-doc-link span {\n",
       "  display: none;\n",
       "  z-index: 9999;\n",
       "  position: relative;\n",
       "  font-weight: normal;\n",
       "  right: .2ex;\n",
       "  padding: .5ex;\n",
       "  margin: .5ex;\n",
       "  width: min-content;\n",
       "  min-width: 20ex;\n",
       "  max-width: 50ex;\n",
       "  color: var(--sklearn-color-text);\n",
       "  box-shadow: 2pt 2pt 4pt #999;\n",
       "  /* unfitted */\n",
       "  background: var(--sklearn-color-unfitted-level-0);\n",
       "  border: .5pt solid var(--sklearn-color-unfitted-level-3);\n",
       "}\n",
       "\n",
       ".sk-estimator-doc-link.fitted span {\n",
       "  /* fitted */\n",
       "  background: var(--sklearn-color-fitted-level-0);\n",
       "  border: var(--sklearn-color-fitted-level-3);\n",
       "}\n",
       "\n",
       ".sk-estimator-doc-link:hover span {\n",
       "  display: block;\n",
       "}\n",
       "\n",
       "/* \"?\"-specific style due to the `<a>` HTML tag */\n",
       "\n",
       "#sk-container-id-1 a.estimator_doc_link {\n",
       "  float: right;\n",
       "  font-size: 1rem;\n",
       "  line-height: 1em;\n",
       "  font-family: monospace;\n",
       "  background-color: var(--sklearn-color-background);\n",
       "  border-radius: 1rem;\n",
       "  height: 1rem;\n",
       "  width: 1rem;\n",
       "  text-decoration: none;\n",
       "  /* unfitted */\n",
       "  color: var(--sklearn-color-unfitted-level-1);\n",
       "  border: var(--sklearn-color-unfitted-level-1) 1pt solid;\n",
       "}\n",
       "\n",
       "#sk-container-id-1 a.estimator_doc_link.fitted {\n",
       "  /* fitted */\n",
       "  border: var(--sklearn-color-fitted-level-1) 1pt solid;\n",
       "  color: var(--sklearn-color-fitted-level-1);\n",
       "}\n",
       "\n",
       "/* On hover */\n",
       "#sk-container-id-1 a.estimator_doc_link:hover {\n",
       "  /* unfitted */\n",
       "  background-color: var(--sklearn-color-unfitted-level-3);\n",
       "  color: var(--sklearn-color-background);\n",
       "  text-decoration: none;\n",
       "}\n",
       "\n",
       "#sk-container-id-1 a.estimator_doc_link.fitted:hover {\n",
       "  /* fitted */\n",
       "  background-color: var(--sklearn-color-fitted-level-3);\n",
       "}\n",
       "</style><div id=\"sk-container-id-1\" class=\"sk-top-container\"><div class=\"sk-text-repr-fallback\"><pre>LogisticRegression()</pre><b>In a Jupyter environment, please rerun this cell to show the HTML representation or trust the notebook. <br />On GitHub, the HTML representation is unable to render, please try loading this page with nbviewer.org.</b></div><div class=\"sk-container\" hidden><div class=\"sk-item\"><div class=\"sk-estimator fitted sk-toggleable\"><input class=\"sk-toggleable__control sk-hidden--visually\" id=\"sk-estimator-id-1\" type=\"checkbox\" checked><label for=\"sk-estimator-id-1\" class=\"sk-toggleable__label fitted sk-toggleable__label-arrow fitted\">&nbsp;&nbsp;LogisticRegression<a class=\"sk-estimator-doc-link fitted\" rel=\"noreferrer\" target=\"_blank\" href=\"https://scikit-learn.org/1.5/modules/generated/sklearn.linear_model.LogisticRegression.html\">?<span>Documentation for LogisticRegression</span></a><span class=\"sk-estimator-doc-link fitted\">i<span>Fitted</span></span></label><div class=\"sk-toggleable__content fitted\"><pre>LogisticRegression()</pre></div> </div></div></div></div>"
      ],
      "text/plain": [
       "LogisticRegression()"
      ]
     },
     "execution_count": 56,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "from sklearn.model_selection import train_test_split\n",
    "from sklearn.linear_model import LogisticRegression\n",
    "\n",
    "X = patients[[\"CAG_Repeats\"]]\n",
    "y = patients[\"Disease\"]\n",
    "\n",
    "X_train, X_test, y_train, y_test = train_test_split(\n",
    "    X, y, test_size=0.2, random_state=42\n",
    ")\n",
    "\n",
    "model = LogisticRegression()\n",
    "model.fit(X_train, y_train)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 58,
   "id": "bcb6a26d-c3c1-425d-9162-54608e20f9fa",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "1.0"
      ]
     },
     "execution_count": 58,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "model.score(X_test, y_test)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 60,
   "id": "59703c89-33b1-4f7e-83e4-efe68c5ee26a",
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "C:\\ProgramData\\anaconda3\\Lib\\site-packages\\sklearn\\base.py:493: UserWarning: X does not have valid feature names, but LogisticRegression was fitted with feature names\n",
      "  warnings.warn(\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "array([1], dtype=int64)"
      ]
     },
     "execution_count": 60,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "model.predict([[40]])"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 64,
   "id": "05c25789-cd1b-479d-bc64-4f8ff1899aa8",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "array([1], dtype=int64)"
      ]
     },
     "execution_count": 64,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "model.predict(pd.DataFrame({\"CAG_Repeats\":[40]}))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 66,
   "id": "8ed1a091-ac3c-41fa-9e14-fddbda2c0ab9",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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XBn4sK+vsm731q1bEEKYJt+WZG1ZaTIyMkp93Rqfv19c974eyyonJwcZGRmm/VpTp04dbNq0CT4+Phg0aBDq1KmDOnXq4OuvvzaVuXbtGm7fvg0HBweL9620tDTcuHEDQPEcsRYtWmD//v2YPHkytm3bhsTERKxcudLseGx5nzK+Vz777LMW+4yPjzftk4pxzo2MxMbGYty4cZg7dy7+/e9/V0oMa9asAQDTREPjHJFvv/221KucfH19odfroVKpLCbO2iIoKAjz588HAJw9exY//fQTJkyYgIKCAsydO9fqOl5eXti/fz+EEGYfhunp6SgqKirz3BZ7GL+tHjt2DJ06dTK94aemplokeVevXjWLSafTWUySBYrvuVKW2I37TktLs3guLS3NbLLqjz/+iODgYMTHx5udO2vxACj1HiN169ZFdHQ0Zs2ahejoaKxZswYTJ04s85wLe+u0Iu99otPpABSfk5JzxR63D6hq1apBoVBYzK8BrL8WSuPl5YXU1FSL5VevXgWAR9KegOI5Rnq9/oH3ZWrRogVatGgBvV6PgwcP4ttvv8WwYcPg6+uLPn36wNvbG15eXli/fr3V9Y29YVu2bMHVq1exbds2U28NAKvz/R70PmU8Rz///LPV3lIyx54bGQkMDMQ///lPdO3aFTExMY98/8eOHcOUKVNQu3Zt9OrVCwDQvHlzeHh44PTp01a/9UZERMDBwQGOjo5o1aoVli9f/lAfAPXr18cnn3yCRo0a4fDhw6WWa9euHe7cuWNxQ7LFixebnq9oxqskjBNZ27ZtC6A4eSgpMTERSUlJZjHVrl0bx48fNyt39uzZMl9x8txzz0Gn05l6UYz27NljMVSgUCjg4OBglhykpaVZvVrqQYYOHYrjx48jJiYGKpUK7777bpniBx6POjUyJoP31tHatWvLvM3y6t0oydnZGREREVi9ejUKCgpMy+/cuWP1qqrSepPatWtn+qAvafHixXBycirT7RvslZycjA8//BDu7u6mCyseRKVSoVmzZpg1axYAmN4zXnzxRWRkZECv11t9zzJOsDe2gXsvdih5tZY11t6nOnXqBLVajfPnz5f6Xkn/w54bmfnss89sLnv+/Hn8/PPPFssbNmz4wFn5hw4dgru7OwoLC0038VuyZAl8fHywdu1a0/i9i4sLvv32W8TExODmzZvo2bMnfHx8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",
      "text/plain": [
       "<Figure size 640x480 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "\n",
    "x_range = np.linspace(15,55,200)\n",
    "pred = model.predict(pd.DataFrame({\"CAG_Repeats\":x_range}))\n",
    "\n",
    "plt.scatter(patients[\"CAG_Repeats\"], patients[\"Disease\"], alpha=0.5)\n",
    "plt.plot(x_range, pred, color=\"red\")\n",
    "plt.xlabel(\"CAG Repeat\")\n",
    "plt.ylabel(\"Disease Risk\")\n",
    "plt.title(\"ML Decision Boundary for Huntington Disease\")\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 68,
   "id": "906a16a3-f58e-44f6-9c22-8ccc5b2ce4b7",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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3qpUqZJuueUWvfrqq/keg9wuG8iud9++fRowYID8/PxUvXp1DR06VKdPny7wuJZFFy5c0IwZM9SkSRONHz8+1zbu7u4aOXJkjvnLly9XmzZtVL58eVWoUEFdu3bNMcazL5X48ccf1b17d1WoUEHBwcF6/PHHlZGRkaOWyZMnO3721apV0wMPPKDffvut+HYYZRLhFZbVvXt3ubm56bPPPsuzzeHDh3XnnXfK09NTb731ltasWaNp06apfPnyunDhggIDA7VmzRpJ0rBhw7Rt2zZt27ZNzzzzjFM/ffr0UZ06dfTee+9p7ty5+da1a9cuPfroo3rssce0atUqhYWF6ZFHHtFLL7101fv4xhtvqG3btgoICHDUtm3btjzbHzhwQGFhYdq3b59ee+01rVy5Uo0aNdKQIUM0Y8aMHO2feuopHTlyRG+++ab+85//6ODBg+rZs6cuXbqUb12bN29Wx44ddfr0aS1YsEDLli2Tr6+vevbsqeXLl0v68yPelStXSvrzY+Zt27Zp1apVefa5bt06ubm5qWfPnoU5NJL+DD3r1q1Tr169VK1aNUVHR+vHH390GhOXLl1SfHy8WrZsqcDAwEL3fTUSExMlSfXq1XPM++6773Tbbbdp7969mjlzpv73v//pzjvv1JgxYxwf1V/uqaee0qFDh/Tmm2/qzTff1C+//KKIiAgdOnTI0SY+Pl5t27bVqVOnNHfuXH344Ydq1qyZ+vXr5xSafv/9d0nSpEmTFBcXp9jYWN10002KiIhw/GFU2LFf1P1PSUnR7bffrrVr1+pf//qXPvnkEw0bNkxTp07VQw895Gh39uxZRUZGKj4+XtOnT9e7776r6tWrq1+/foXa9sMPP6wnn3xSnTt31urVq/X8889rzZo1CgsL0/Hjx53apqSk6L777tP999+v1atXKyoqShMmTHBcdiJJM2bMUExMjAYMGKC4uDgtX75cw4YN06lTp676uGS75557VK9ePa1YsULjx4/X0qVL9dhjjxV6/YsXL+Y6ucKOHTt06tSpq/o9lf78A3LAgAFq1KiR3n33Xf33v/9Venq62rdvr++++86pbWZmpu666y516tRJH374oYYOHaqXX35Z06dPd7TJyspSr169NG3aNA0cOFBxcXGaNm2a1q9fr4iICP3xxx/Fsr8oo1x96hfIS2E+zqpevbpp2LCh4/WkSZPM5cP6/fffN5LMrl278uwjv49Os/vL7ePgK7dlzJ+XDdhsthzb69y5s6lYsaLjkoPsfUtMTHRql/2xdHx8vGNefpcNXFl3//79jd1uN0lJSU7toqKijI+Pjzl16pTTdrp37+7U7t133zWSzLZt23LdXrbWrVsbf39/k56e7ph38eJF06RJE1OjRg2TlZVljMn7Y+bcNGjQwAQEBBTY7nLPPfeckWTWrFljjDHm0KFDxmazmUGDBjnapKSkGEmmf//+Oda/ePGiyczMdEzZdecl+7gtX77cZGZmmnPnzpkvvvjC1K9f3zRq1MicPHnS0bZr166mRo0a5vTp0059jBo1ynh5eZnff//dqc8WLVo4bf/w4cPGw8PDPPjgg455DRo0MM2bNzeZmZlOffbo0cMEBgaaS5cu5Vp39n526tTJ9O7d2zG/qJcNTJ8+3WRmZprz58+bXbt2mTZt2pjAwECn8Tx8+HBToUIFc+TIEac+XnrpJSPJ7Nu3zxhjzL///W8jyXzyySdO7YYPH24kmdjYWMe8K3/n9u/fbySZESNGOK375ZdfGknmqaeecswLDw83ksyXX37p1LZRo0ama9eujtc9evQwzZo1K9TxuFz2scmt3hkzZji1HTFihPHy8ipwvEVHRxtJ+U6lfdnAO++8YySZuXPn5lh2+e/S5WM0KSnJuLu7m9GjRzu1T09PNwEBAaZv37459vndd991atu9e3dTv359x+tly5YZSWbFihVO7RISEowk88YbbxRqf2BNnHmFpZkrPh6+UrNmzeTp6am///3vWrhwodNZrKtxzz33FLpt48aN1bRpU6d5AwcOVFpamr755psibb+wNm7cqE6dOik4ONhp/pAhQ3Tu3LkcZ23vuusup9e33HKLJOnIkSN5buPs2bP68ssvde+99zrdCe3m5qZBgwbp2LFjhb704K8wxjguFejcubMkKTQ0VBEREVqxYoXTzSt5ufXWW+Xh4eGYZs6cWaht9+vXTx4eHvLx8VHbtm2VlpamuLg4VapUSZJ0/vx5ffrpp+rdu7d8fHyczpZ1795d58+f1/bt2536HDhwoNNH4iEhIQoLC1N8fLwk6ccff9T333+v++67T5Jy9JmcnOx03OfOnasWLVrIy8tL7u7u8vDw0Keffqr9+/cXah/z8+STT8rDw8NxWcrevXv10UcfOV0W8r///U+RkZEKCgpyqjUqKkrSn2fvs//r6+ub40asAQMGFFhH9rG58lKH22+/XQ0bNtSnn37qND8gIEC3336707xbbrnFabzffvvt2r17t0aMGKG1a9cWahwVJLffs/Pnz+f6xJQreXt7KyEhIdfJ29v7L9dWXE6dOuX0u+Th4eF4AsfatWt18eJFDR482GkseHl5KTw8PMfTLmw2W44zu1f+nP73v/+pUqVK6tmzp1OfzZo1U0BAwFU/PQbWQniFZZ09e1YnTpzI9zrF2rVra8OGDfL399fIkSNVu3Zt1a5du8Dr1650NR83BwQE5DnvxIkTV7Xdq3XixIlca80+Rlduv0qVKk6v7Xa7JOX7kdvJkydljLmq7RRGzZo19dtvv+ns2bOFar9x40YlJibqb3/7m9LS0nTq1CmdOnVKffv21blz57Rs2TJJUtWqVeXt7Z1rIF+6dKkSEhK0evXqq6p1+vTpSkhI0ObNmzVx4kT9+uuvuvvuux3X5J04cUIXL17U66+/nuMf9O7du0tSjo+08xo32ccy+3rScePG5ehzxIgRTn3OmjVLDz/8sFq1aqUVK1Zo+/btSkhIULdu3Yrl49RHHnlECQkJ2rJli1566SVlZmaqV69eTj/3X3/9VR999FGOWhs3buxU64kTJ1S9evUc28ht3pWyt5fXWCxovEt/jvnLj8mECRP00ksvafv27YqKilKVKlXUqVOnv/QotKL8nmUrV66cWrZsmetUrlzp/xNes2ZNSTn/wPX19XWE6kmTJjktyx67t912W47xsHz58hy/Cz4+PvLy8nKaZ7fbdf78eac+T506JU9Pzxx9pqSk5OgT1xYelQXLiouL06VLlwp8xE/79u3Vvn17Xbp0STt27NDrr7+uRx99VNWrV1f//v0Lta2rebZkSkpKnvOy/xHLfmO+8gaEv/qGW6VKFSUnJ+eY/8svv0j6M8j9VZUrV1a5cuWKfTtdu3bVunXr9NFHHxXq57JgwQJJfwa1WbNm5bp8+PDhcnNzU8eOHbVu3TolJyc7BZ1GjRpJUo4b5wpy0003OW5S6tChg7y9vfX000/r9ddf17hx41S5cmXHmejcblyR/jxLfLm8xk32mMk+phMmTFCfPn1y7TP7MWSLFy9WRESE5syZ47Q8PT39KvYybzVq1HDsf/Y12ffff78mTZqk2bNnO+q95ZZb9MILL+TaR/YfOlWqVNFXX32VY3lux+NK2ccmOTk5x01+v/zyS5HGobu7u8aOHauxY8fq1KlT2rBhg5566il17dpVR48elY+Pz1X3eS259dZbVblyZX300UeaMmWKY76bm5tjTOzdu9dpneyfw/vvv6+QkJBiqaNq1aqqUqWK47rtK/n6+hbLdlA2ceYVlpSUlKRx48bJz89Pw4cPL9Q6bm5uatWqleNu7+yP8K/mLEhh7Nu3T7t373aat3TpUvn6+qpFixaS5Ph4dc+ePU7tcjsDeOWZofx06tRJGzdudITIbIsWLZKPj49at25d2N3IU/ny5dWqVSutXLnSqa6srCwtXrxYNWrUcLpxp7CGDRumgIAAPfHEE/r5559zbZN9A9jJkye1atUqtW3bVvHx8Tmm7LvVs/8RnTBhgi5duqR//OMfyszMLMJe5++JJ55QnTp1NG3aNKWnp8vHx0eRkZHauXOnbrnlllzPml15Nm7ZsmVOl8EcOXJEW7dudfxxVr9+fdWtW1e7d+/O80xc9j/YNpvNMa6z7dmzJ8dlI8U19u+77z5FRERo/vz5jjNyPXr00N69e1W7du1ca80Or+Hh4UpPT9cnn3zi1Oc777xT4HY7duwoSU43XElSQkKC9u/fr06dOv2l/apUqZLuvfdejRw5Ur///vtV/5FzLfL09NQ///lP7d271+kGqvx07dpV7u7u+umnn/Icu1erR48eOnHihC5dupRrf/k9TxrWx5lXlHl79+51XM+Umpqqzz//XLGxsXJzc9OqVascj/vJzdy5c7Vx40bdeeedqlmzps6fP6+33npLkhxfbuDr66uQkBB9+OGH6tSpk2644QZVrVo138c65ScoKEh33XWXYmJiFBgYqMWLF2v9+vWaPn2646zNbbfdpvr162vcuHG6ePGiKleurFWrVmnLli05+rv55pu1cuVKzZkzR7feeqvjY8TcTJo0yXGt4b/+9S/dcMMNWrJkieLi4jRjxgz5+fkVaZ+uNHXqVHXu3FmRkZEaN26cPD099cYbb2jv3r1atmxZkb4Fyc/PTx9++KF69Oih5s2bO31JwcGDB7V48WLt3r1bffr00ZIlS3T+/HmNGTMm1zPvVapU0ZIlS7RgwQK9/PLLatu2rf79739r9OjRatGihf7+97+rcePGjjPIK1askCRVrFixSMfDw8NDU6ZMUd++ffXqq6/q6aef1quvvqp27dqpffv2evjhh1WrVi2lp6frxx9/1EcffaSNGzc69ZGamqrevXvroYce0unTpzVp0iR5eXlpwoQJjjbz5s1TVFSUunbtqiFDhujGG2/U77//rv379+ubb77Re++9J+nPf9iff/55TZo0SeHh4Tpw4ICee+45hYaGOt2lXpxjf/r06WrVqpWef/55vfnmm3ruuee0fv16hYWFacyYMapfv77Onz+vw4cP6+OPP9bcuXNVo0YNRUdH6+WXX9b999+vyZMnq06dOvrkk0+0du1aScr3o/H69evr73//u15//XWVK1dOUVFROnz4sJ555hkFBwdf1R392Xr27KkmTZqoZcuWqlatmo4cOaJXXnlFISEhlvoCipL05JNP6vvvv9f48eP12WefqV+/fqpVq5YyMjIcT8xwc3NzvN/VqlVLzz33nCZOnKhDhw6pW7duqly5sn799Vd99dVXKl++fK5P4MhP//79tWTJEnXv3l2PPPKIbr/9dnl4eOjYsWOKj49Xr1691Lt375LYfZQFrr1fDMhb9l2w2ZOnp6fx9/c34eHhZsqUKSY1NTXHOlfejbxt2zbTu3dvExISYux2u6lSpYoJDw83q1evdlpvw4YNpnnz5sZutzvdwZvd32+//Vbgtoz5/19S8P7775vGjRsbT09PU6tWLTNr1qwc6//www+mS5cupmLFiqZatWpm9OjRJi4uLsfTBn7//Xdz7733mkqVKhmbzea0TeVyp/i3335revbsafz8/Iynp6dp2rSp0x3QxuT9sP3c7pjOy+eff246duxoypcvb7y9vU3r1q3NRx99lGt/hXnaQLaUlBTz5JNPmsaNGxsfHx9jt9tNnTp1zPDhw823335rjDGmWbNmxt/fP88vHDDmzyciVK1a1anNrl27zAMPPGBCQ0ON3W43Xl5epk6dOmbw4MHm008/LbC2vI5btlatWpnKlSs7nuqQmJhohg4dam688Ubj4eFhqlWrZsLCwszkyZNz9Pnf//7XjBkzxlSrVs3Y7XbTvn17s2PHjhzb2L17t+nbt6/x9/c3Hh4eJiAgwHTs2NHp7u+MjAwzbtw4c+ONNxovLy/TokUL88EHH5jo6OgcT67Ia+znpqCf59/+9jfj7u5ufvzxR2PMn08zGDNmjAkNDTUeHh7mhhtuMLfeequZOHGiOXPm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NWuWfv75Z8XGxiolJUXt2rVT48aNNWvWLP3222/FXScAAADw127Ycnd3V+/evfXuu+9q+vTp+umnnzRu3DjVqFFDgwcPVnJycnHVCQAAAPy18Lpjxw6NGDFCgYGBmjVrlsaNG6effvpJGzdu1M8//6xevXoVV50AAACA3Iuy0qxZsxQbG6sDBw6oe/fuWrRokbp3765y5f7MwqGhoZo3b54aNGhQrMUCAADg+lak8DpnzhwNHTpUDzzwgAICAnJtU7NmTS1YsOAvFQcAAABcrkjh9eDBgwW28fT0VHR0dFG6BwAAAHJVpGteY2Nj9d577+WY/95772nhwoV/uSgAAAAgN0UKr9OmTVPVqlVzzPf399eUKVP+clEAAABAbooUXo8cOaLQ0NAc80NCQpSUlPSXiwIAAAByU6Tw6u/vrz179uSYv3v3blWpUuUvFwUAAADkpkjhtX///hozZozi4+N16dIlXbp0SRs3btQjjzyi/v37F3eNAAAAgKQiPm1g8uTJOnLkiDp16iR39z+7yMrK0uDBg7nmFQAAACWmSOHV09NTy5cv1/PPP6/du3fL29tbN998s0JCQoq7PgAAAMChSOE1W7169VSvXr3iqgUAAADIV5Gueb106ZIWLFiggQMH6o477lDHjh2dpsL67LPP1LNnTwUFBclms+mDDz5wWj5kyBDZbDanqXXr1kUpGQAAANeAIp15feSRR/T222/rzjvvVJMmTWSz2Yq08bNnz6pp06Z64IEHdM899+Taplu3boqNjXW89vT0LNK2AAAAYH1FCq/vvPOO3n33XXXv3v0vbTwqKkpRUVH5trHb7QoICPhL2wEAAMC1oUiXDXh6eqpOnTrFXUuuNm3aJH9/f9WrV08PPfSQUlNT822fkZGhtLQ0pwkAAADXhiKF18cff1yvvvqqjDHFXY+TqKgoLVmyRBs3btTMmTOVkJCgjh07KiMjI891pk6dKj8/P8cUHBxcojUCAACg9BTpsoEtW7YoPj5en3zyiRo3biwPDw+n5StXriyW4vr16+f4/yZNmqhly5YKCQlRXFyc+vTpk+s6EyZM0NixYx2v09LSCLAAAADXiCKF10qVKql3797FXUuBAgMDFRISooMHD+bZxm63y263l2JVAAAAKC1FCq+X3/1fmk6cOKGjR48qMDDQJdsHAACAaxXpmldJunjxojZs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kv48AxsYGOiYn5qamuNs7OXsdrvsdnuJ1wcAAIDSV2bPvIaGhiogIEDr1693zLtw4YI2b96ssLAwF1YGAAAAV3HpmdczZ87oxx9/dLxOTEzUrl27dMMNN6hmzZp69NFHNWXKFNWtW1d169bVlClT5OPjo4EDB7qwagAAALiKS8Prjh07FBkZ6Xg9duxYSVJ0dLTefvttPfHEE/rjjz80YsQInTx5Uq1atdK6devk6+vrqpIBAADgQi4NrxERETLG5LncZrMpJiZGMTExpVcUAAAAyqwye80rAAAAcCXCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsIwyHV5jYmJks9mcpoCAAFeXBQAAABdxd3UBBWncuLE2bNjgeO3m5ubCagAAAOBKZT68uru7c7YVAAAAksr4ZQOSdPDgQQUFBSk0NFT9+/fXoUOH8m2fkZGhtLQ0pwkAAADXhjIdXlu1aqVFixZp7dq1mj9/vlJSUhQWFqYTJ07kuc7UqVPl5+fnmIKDg0uxYgAAAJSkMh1eo6KidM899+jmm2/WHXfcobi4OEnSwoUL81xnwoQJOn36tGM6evRoaZULAACAElbmr3m9XPny5XXzzTfr4MGDebax2+2y2+2lWBUAAABKS5k+83qljIwM7d+/X4GBga4uBQAAAC5QpsPruHHjtHnzZiUmJurLL7/Uvffeq7S0NEVHR7u6NAAAALhAmb5s4NixYxowYICOHz+uatWqqXXr1tq+fbtCQkJcXRoAAABcoEyH13feecfVJQAAAKAMKdOXDQAAAACXI7wCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMgivAAAAsAzCKwAAACyD8AoAAADLILwCAADAMtxdXQAAAIDVJCUl6fjx464uo0RVrVpVNWvWdHUZORBeAQAArkJSUpIaNGyoP86dc3UpJcrbx0ff799f5gIs4RUAAOAqHD9+XH+cO6e+k+fIP7Suq8spEamJB/Xu0w/r+PHjhFcAAIBrgX9oXd3YsKmry7jucMMWAAAALIPwCgAAAMsgvAIAAMAyCK8AAACwDMIrAAAALIPwCgAAAMsgvAIAAMAyLBFe33jjDYWGhsrLy0u33nqrPv/8c1eXBAAAABco8+F1+fLlevTRRzVx4kTt3LlT7du3V1RUlJKSklxdGgAAAEpZmQ+vs2bN0rBhw/Tggw+qYcOGeuWVVxQcHKw5c+a4ujQAAACUsjL99bAXLlzQ119/rfHjxzvN79Kli7Zu3ZrrOhkZGcrIyHC8Pn36tCQpLS2t5Aq9zJkzZyRJP+/fowvnzpbKNkvbb0d+kiR9/fXXjv291pQrV05ZWVmuLqPEXOv7J137+8j+Wd+1vo/X8v4dOHBA0vXxb/2ZM2dKJUNlb8MYU3BjU4b9/PPPRpL54osvnOa/8MILpl69ermuM2nSJCOJiYmJiYmJiYnJYtPRo0cLzIdl+sxrNpvN5vTaGJNjXrYJEyZo7NixjtdZWVn6/fffVaVKlTzXKSvS0tIUHByso0ePqmLFiq4up8ziOBUOx6nwOFaFw3EqHI5T4XCcCud6OU7GGKWnpysoKKjAtmU6vFatWlVubm5KSUlxmp+amqrq1avnuo7dbpfdbneaV6lSpZIqsURUrFjxmh6gxYXjVDgcp8LjWBUOx6lwOE6Fw3EqnOvhOPn5+RWqXZm+YcvT01O33nqr1q9f7zR//fr1CgsLc1FVAAAAcJUyfeZVksaOHatBgwapZcuWatOmjf7zn/8oKSlJ//jHP1xdGgAAAEpZmQ+v/fr104kTJ/Tcc88pOTlZTZo00ccff6yQkBBXl1bs7Ha7Jk2alOOyBzjjOBUOx6nwOFaFw3EqHI5T4XCcCofjlJPNmMI8kwAAAABwvTJ9zSsAAABwOcIrAAAALIPwCgAAAMsgvAIAAMAyCK+lZOrUqbrtttvk6+srf39/3X333Y7vRs7Lpk2bZLPZckzff/99KVVd+mJiYnLsb0BAQL7rbN68Wbfeequ8vLx00003ae7cuaVUrevUqlUr17ExcuTIXNtfT2Pps88+U8+ePRUUFCSbzaYPPvjAabkxRjExMQoKCpK3t7ciIiK0b9++AvtdsWKFGjVqJLvdrkaNGmnVqlUltAelI7/jlJmZqSeffFI333yzypcvr6CgIA0ePFi//PJLvn2+/fbbuY6z8+fPl/DelJyCxtOQIUNy7G/r1q0L7Pd6Gk+Sch0XNptNL774Yp59XovjqTBZgPeoghFeS8nmzZs1cuRIbd++XevXr9fFixfVpUsXnT17tsB1Dxw4oOTkZMdUt27dUqjYdRo3buy0v99++22ebRMTE9W9e3e1b99eO3fu1FNPPaUxY8ZoxYoVpVhx6UtISHA6Rtlf5PG3v/0t3/Wuh7F09uxZNW3aVLNnz851+YwZMzRr1izNnj1bCQkJCggIUOfOnZWenp5nn9u2bVO/fv00aNAg7d69W4MGDVLfvn315ZdfltRulLj8jtO5c+f0zTff6JlnntE333yjlStX6ocfftBdd91VYL8VK1Z0GmPJycny8vIqiV0oFQWNJ0nq1q2b0/5+/PHH+fZ5vY0nSTnGxFtvvSWbzaZ77rkn336vtfFUmCzAe1QhGLhEamqqkWQ2b96cZ5v4+HgjyZw8ebL0CnOxSZMmmaZNmxa6/RNPPGEaNGjgNG/48OGmdevWxVxZ2fbII4+Y2rVrm6ysrFyXX49jyRhjJJlVq1Y5XmdlZZmAgAAzbdo0x7zz588bPz8/M3fu3Dz76du3r+nWrZvTvK5du5r+/fsXe82ucOVxys1XX31lJJkjR47k2SY2Ntb4+fkVb3FlSG7HKTo62vTq1euq+mE8GdOrVy/TsWPHfNtc6+PJmJxZgPeowuHMq4ucPn1aknTDDTcU2LZ58+YKDAxUp06dFB8fX9KludzBgwcVFBSk0NBQ9e/fX4cOHcqz7bZt29SlSxeneV27dtWOHTuUmZlZ0qWWCRcuXNDixYs1dOhQ2Wy2fNteb2PpSomJiUpJSXEaM3a7XeHh4dq6dWue6+U1zvJb51pz+vRp2Ww2VapUKd92Z86cUUhIiGrUqKEePXpo586dpVOgC23atEn+/v6qV6+eHnroIaWmpubb/nofT7/++qvi4uI0bNiwAtte6+PpyizAe1ThEF5dwBijsWPHql27dmrSpEme7QIDA/Wf//xHK1as0MqVK1W/fn116tRJn332WSlWW7patWqlRYsWae3atZo/f75SUlIUFhamEydO5No+JSVF1atXd5pXvXp1Xbx4UcePHy+Nkl3ugw8+0KlTpzRkyJA821yPYyk3KSkpkpTrmMleltd6V7vOteT8+fMaP368Bg4cqIoVK+bZrkGDBnr77be1evVqLVu2TF5eXmrbtq0OHjxYitWWrqioKC1ZskQbN27UzJkzlZCQoI4dOyojIyPPda738bRw4UL5+vqqT58++ba71sdTblmA96jCKfNfD3stGjVqlPbs2aMtW7bk265+/fqqX7++43WbNm109OhRvfTSS+rQoUNJl+kSUVFRjv+/+eab1aZNG9WuXVsLFy7U2LFjc13nyrON5v++NK6gs5DXigULFigqKkpBQUF5trkex1J+chszBY2XoqxzLcjMzFT//v2VlZWlN954I9+2rVu3drpZqW3btmrRooVef/11vfbaayVdqkv069fP8f9NmjRRy5YtFRISori4uHzD2fU6niTprbfe0n333VfgtavX+njKLwvwHpU/zryWstGjR2v16tWKj49XjRo1rnr91q1bXzN/dRZG+fLldfPNN+e5zwEBATn+skxNTZW7u7uqVKlSGiW61JEjR7RhwwY9+OCDV73u9TaWJDmeXJHbmLnyrMWV613tOteCzMxM9e3bV4mJiVq/fn2+Z11zU65cOd12223X1TgLDAxUSEhIvvt8vY4nSfr888914MCBIr1nXUvjKa8swHtU4RBeS4kxRqNGjdLKlSu1ceNGhYaGFqmfnTt3KjAwsJirK7syMjK0f//+PPe5TZs2jjvts61bt04tW7aUh4dHaZToUrGxsfL399edd9551eteb2NJkkJDQxUQEOA0Zi5cuKDNmzcrLCwsz/XyGmf5rWN12cH14MGD2rBhQ5H+GDTGaNeuXdfVODtx4oSOHj2a7z5fj+Mp24IFC3TrrbeqadOmV73utTCeCsoCvEcVkmvuE7v+PPzww8bPz89s2rTJJCcnO6Zz58452owfP94MGjTI8frll182q1atMj/88IPZu3evGT9+vJFkVqxY4YpdKBWPP/642bRpkzl06JDZvn276dGjh/H19TWHDx82xuQ8RocOHTI+Pj7mscceM999951ZsGCB8fDwMO+//76rdqHUXLp0ydSsWdM8+eSTOZZdz2MpPT3d7Ny50+zcudNIMrNmzTI7d+503CU/bdo04+fnZ1auXGm+/fZbM2DAABMYGGjS0tIcfQwaNMiMHz/e8fqLL74wbm5uZtq0aWb//v1m2rRpxt3d3Wzfvr3U96+45HecMjMzzV133WVq1Khhdu3a5fSelZGR4ejjyuMUExNj1qxZY3766Sezc+dO88ADDxh3d3fz5ZdfumIXi0V+xyk9Pd08/vjjZuvWrSYxMdHEx8ebNm3amBtvvJHxdMXvnTHGnD592vj4+Jg5c+bk2sf1MJ4KkwV4jyoY4bWUSMp1io2NdbSJjo424eHhjtfTp083tWvXNl5eXqZy5cqmXbt2Ji4urvSLL0X9+vUzgYGBxsPDwwQFBZk+ffqYffv2OZZfeYyMMWbTpk2mefPmxtPT09SqVSvPN8Zrzdq1a40kc+DAgRzLruexlP1YsCun6OhoY8yfj6KZNGmSCQgIMHa73XTo0MF8++23Tn2Eh4c72md77733TP369Y2Hh4dp0KCB5YN/fscpMTExz/es+Ph4Rx9XHqdHH33U1KxZ03h6eppq1aqZLl26mK1bt5b+zhWj/I7TuXPnTJcuXUy1atWMh4eHqVmzpomOjjZJSUlOfVzv4ynbvHnzjLe3tzl16lSufVwP46kwWYD3qILZjPm/u1sAAACAMo5rXgEAAGAZhFcAAABYBuEVAAAAlkF4BQAAgGUQXgEAAGAZhFcAAABYBuEVAAAAlkF4BQAAgGUQXgEAJebw4cOy2WzatWuXq0sBcI0gvAK4rgwZMkQ2m002m03u7u6qWbOmHn74YZ08edLVpeVQ2OBXVgLikCFDdPfdd7u0BgDXPsIrgOtOt27dlJycrMOHD+vNN9/URx99pBEjRri6LABAIRBeAVx37Ha7AgICVKNGDXXp0kX9+vXTunXrnNrExsaqYcOG8vLyUoMGDfTGG284lmWf6XznnXcUFhYmLy8vNW7cWJs2bXLq47vvvlP37t1VoUIFVa9eXYMGDdLx48cdy9esWaN27dqpUqVKqlKlinr06KGffvrJsTw0NFSS1Lx5c9lsNkVERBRpf40xmjFjhm666SZ5e3uradOmev/99x3LN23aJJvNpk8//VQtW7aUj4+PwsLCdODAAad+Jk+eLH9/f/n6+urBBx/U+PHj1axZM0lSTEyMFi5cqA8//NBxZvvy43Ho0CFFRkbKx8dHTZs21bZt24q0LwBAeAVwXTt06JDWrFkjDw8Px7z58+dr4sSJeuGFF7R//35NmTJFzzzzjBYuXOi07j//+U89/vjj2rlzp8LCwnTXXXfpxIkTkqTk5GSFh4erWbNm2rFjh9asWaNff/1Vffv2dax/9uxZjR07VgkJCfr0009Vrlw59e7dW1lZWZKkr776SpK0YcMGJScna+XKlUXax6efflqxsbGaM2eO9u3bp8cee0z333+/Nm/e7NRu4sSJmjlzpnbs2CF3d3cNHTrUsWzJkiV64YUXNH36dH399deqWbOm5syZ41g+btw49e3b13FWOzk5WWFhYU59jxs3Trt27VK9evU0YMAAXbx4sUj7A+A6ZwDgOhIdHW3c3NxM+fLljZeXl5FkJJlZs2Y52gQHB5ulS5c6rff888+bNm3aGGOMSUxMNJLMtGnTHMszMzNNjRo1zPTp040xxjzzzDOmS5cuTn0cPXrUSDIHDhzItbbU1FQjyXz77bdO29m5c2e++5RfuzNnzhgvLy+zdetWp/nDhg0zAwYMMMYYEx8fbySZDRs2OJbHxcUZSeaPP/4wxhjTqlUrM3LkSKc+2rZta5o2bep4HR0dbXr16pVrbW+++aZj3r59+4wks3///nz3CwByw5lXANedyMhI7dq1S19++aVGjx6trl27avTo0ZKk3377TUePHtWwYcNUoUIFxzR58mSnj/QlqU2bNo7/d3d3V8uWLbV//35J0tdff634+HinPho0aCBJjn5++uknDRw4UDfddJMqVqzouEwgKSmp2Pb1u+++0/nz59W5c2enWhYtWpRjf2655RbH/wcGBkqSUlNTJUkHDhzQ7bff7tT+ytf5ya9vALga7q4uAABKW/ny5VWnTh1J0muvvabIyEg9++yzev755x0f2c+fP1+tWrVyWs/Nza3Avm02myQpKytLPXv21PTp03O0yQ5vPXv2VHBwsObPn6+goCBlZWWpSZMmunDhwl/av8tl709cXJxuvPFGp2V2u93p9eWXTly+H1fOy2aMKXQdBfUNAIVFeAVw3Zs0aZKioqL08MMPKygoSDfeeKMOHTqk++67L9/1tm/frg4dOkiSLl68qK+//lqjRo2SJLVo0UIrVqxQrVq15O6e8632xIkT2r9/v+bNm6f27dtLkrZs2eLUxtPTU5J06dKlIu9bo0aNZLfblZSUpPDw8CL3U79+fX311VcaNGiQY96OHTuc2nh6ev6lWgGgMAivAK57ERERaty4saZMmaLZs2crJiZGY8aMUcWKFRUVFaWMjAzt2LFDJ0+e1NixYx3r/fvf/1bdunXVsGFDvfzyyzp58qTjJqeRI0dq/vz5GjBggP75z3+qatWq+vHHH/XOO+9o/vz5qly5sqpUqaL//Oc/CgwMVFJSksaPH+9Ul7+/v7y9vbVmzRrVqFFDXl5e8vPzy3M/rnw6gPRneB03bpwee+wxZWVlqV27dkpLS9PWrVtVoUIFRUdHF+oYjR49Wg899JBatmypsLAwLV++XHv27NFNN93kaFOrVi2tXbtWBw4cUJUqVfKtFQCKzNUX3QJAacrtpiJjjFmyZInx9PQ0SUlJjtfNmjUznp6epnLlyqZDhw5m5cqVxpj/fxPS0qVLTatWrYynp6dp2LCh+fTTT536/OGHH0zv3r1NpUqVjLe3t2nQoIF59NFHTVZWljHGmPXr15uGDRsau91ubrnlFrNp0yYjyaxatcrRx/z5801wcLApV66cCQ8Pz3WfsuvJbUpMTDRZWVnm1VdfNfXr1zceHh6mWrVqpmvXrmbz5s3GmP9/w9bJkycdfe7cudOxfrbnnnvOVK1a1VSoUMEMHTrUjBkzxrRu3dqxPDU11XTu3NlUqFDBSDLx8fG53kx28uRJx3IAuFo2Y67ioiUAgA4fPqzQ0FDt3LnT8ZzT61Hnzp0VEBCg//73v64uBcB1hMsGAAAFOnfunObOnauuXbvKzc1Ny5Yt04YNG7R+/XpXlwbgOkN4BQAUyGaz6eOPP9bkyZOVkZGh+vXra8WKFbrjjjtcXRqA6wyXDQAAAMAy+JICAAAAWAbhFQAAAJZBeAUAAIBlEF4BAABgGYRXAAAAWAbhFQAAAJZBeAUAAIBlEF4BAABgGf8Pd3b0uY4ELVoAAAAASUVORK5CYII=",
      "text/plain": [
       "<Figure size 800x500 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.figure(figsize=(8,5))\n",
    "plt.hist(repeat_counts, bins=10, color=\"skyblue\", edgecolor=\"black\")\n",
    "plt.title(\"Distribution of CAG Repeat Regions in HTT Gene\")\n",
    "plt.xlabel(\"Repeat Length\")\n",
    "plt.ylabel(\"Frequency\")\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 70,
   "id": "b420a967-6381-4e1d-b640-9915dfa2307c",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAA1cAAAGHCAYAAABcY6j2AAAAOXRFWHRTb2Z0d2FyZQBNYXRwbG90bGliIHZlcnNpb24zLjkuMiwgaHR0cHM6Ly9tYXRwbG90bGliLm9yZy8hTgPZAAAACXBIWXMAAA9hAAAPYQGoP6dpAABTlElEQVR4nO3dd3wVVf7/8fcl5SaE5NJJAoGAIB2kCQSkSBdRFhEEpahro4MIgrrg6hLB1R+CgosiiKjASpEVpSlFligdVIqUQFAT6Qk1kOT8/uCbu1xSb5gkN+H1fDzmIffMmTOfM2dyzSczc8ZmjDECAAAAANySIvkdAAAAAAAUBiRXAAAAAGABkisAAAAAsADJFQAAAABYgOQKAAAAACxAcgUAAAAAFiC5AgAAAAALkFwBAAAAgAVIrgAAAADAAiRXADzenj179Pjjj6ty5cry8/NTsWLF1LBhQ02ZMkVnzpxJd5uGDRvKZrPpn//8Z6Ztf/XVV3rwwQcVGhoqX19fBQYGqkGDBpowYYJiYmKyjG3ixImy2WzOxcfHRxUrVtRTTz2luLi4HPU3P+3du1cTJ07U0aNHs1V/7ty5stls2rZtW7rr77//foWHh+colhkzZmju3Lk52jY3fPbZZ5o6dWp+h1EgrV+/XjabTV988UW664cMGSKbzSbpf+dUVkt4eHi262UlOjpaw4YNU82aNRUQECA/Pz+Fh4frscce07p162SMsfJwACjEvPM7AADIzAcffKBBgwapevXqeuGFF1SrVi1du3ZN27Zt0/vvv6+oqCgtXbrUZZtdu3Zp586dkqTZs2dr9OjRadpNSUnR448/rnnz5qlLly6KjIxUeHi4Ll++rK1bt2rOnDn66KOPdPz48WzFuXLlSjkcDl24cEGrV6/WW2+9pc2bN2vXrl3y8fG59QORR/bu3atXX31Vbdq0yXFSZJUZM2aodOnSGjhwYL7Gkeqzzz7Tzz//rBEjRuR3KIVa165dFRUV5VLWvHlz9ezZU88//7yzzMvLS8nJyVnWs9vtme5v+fLl6tu3r0qXLq1nn31WDRs2lN1u16FDh/TFF1/o3nvv1dq1a9WuXTsLegegsCO5AuCxoqKi9Nxzz6lDhw5atmyZyy9JHTp00PPPP6+VK1em2e7DDz+UdP2XtBUrVmjz5s2KiIhwqTN58mTNmzdPkZGRevHFF13Wde7cWePGjdO//vWvbMfaqFEjlS5dWpLUvn17nTp1SnPmzNGmTZvUtm3bbLcDWOnSpUsqWrRofofhljJlyqhMmTJpysuVK6dmzZpluX1260nS4cOH1adPH9WuXVtr165VUFCQc13r1q315JNPav369SpRokT2OwDgtsZtgQA81qRJk2Sz2TRr1qx0//rs6+urBx54wKXsypUr+uyzz9SoUSP9v//3/yRJH330kUudq1evasqUKapTp06axCqVt7e3Bg8enOPYGzduLEn6888/XcpT/wIeFBSkokWLqkWLFvr2229d6qTearhz50716NFDQUFBcjgceuyxx3Ty5Mk0+1q4cKGaN2+ugIAAFStWTJ06dXJeuUu1bds2PfLIIwoPD5e/v7/Cw8PVp08fHTt2zFln7ty5evjhhyVJbdu2dd5WZfWteVeuXNG4ceNUuXJl+fr6qnz58ho8eLDOnTvnrBMeHq5ffvlFGzZsSHN7V0pKil5//XVVr15d/v7+Kl68uOrVq6d33nnHub07xzAlJUVTpkxRjRo1ZLfbVbZsWfXv31+//fabs06bNm20YsUKHTt2zOWWs1QzZ85U/fr1VaxYMQUGBqpGjRoaP358lsfi1VdfVdOmTVWyZEkFBQWpYcOGmj17drq3oX322Wdq3ry5ihUrpmLFiumuu+7S7NmzXWKsU6eONm7cqIiICBUtWlRPPPGEJCkmJkaPPfaYypYtK7vdrpo1a+qtt95SSkqKyz6y6selS5c0evRo5y26JUuWVOPGjfX5559n2VdP9Pbbb+vSpUuaMWOGS2J1ozZt2qh+/fouZQcPHlTfvn1djud7773nUif1VsjPP/9cL730kkJDQxUUFKT27dvrwIEDafaTne8GAJ6P5AqAR0pOTtZ3332nRo0aKSwsLNvbLVmyRGfPntUTTzyhatWqqWXLllq4cKEuXLjgrLNt2zadO3dO3bp1y43QJV1/hkOS7rzzTmfZ/Pnz1bFjRwUFBenjjz/WokWLVLJkSXXq1CndX6L+8pe/qGrVqvriiy80ceJELVu2TJ06ddK1a9ecdSZNmqQ+ffqoVq1aWrRokT755BOdP39e99xzj/bu3eusd/ToUVWvXl1Tp07VqlWrNHnyZMXGxqpJkyY6deqUpOtX+iZNmiRJeu+99xQVFaWoqCh17do1y/4mJycrKSkpzXJzkmCMUffu3fXPf/5T/fr104oVKzRq1Ch9/PHHuvfee5WYmChJWrp0qapUqaIGDRo440i9/XPKlCmaOHGi+vTpoxUrVmjhwoV68sknXZIzd47hc889p7Fjx6pDhw5avny5XnvtNa1cuVIRERHOYzNjxgy1aNFCwcHBznhSb11bsGCBBg0apNatW2vp0qVatmyZRo4cqYsXL2Z53I4ePapnnnlGixYt0pIlS9SjRw8NHTpUr732mku9v/3tb3r00UcVGhqquXPnaunSpRowYIBLcixJsbGxeuyxx9S3b199/fXXGjRokE6ePKmIiAitXr1ar732mpYvX6727dtr9OjRGjJkiHPb7PRj1KhRmjlzpoYNG6aVK1fqk08+0cMPP6zTp09n2VfpeiKbnfMkr6xZs0YhISHOP4Zkx969e9WkSRP9/PPPeuutt/TVV1+pa9euGjZsmF599dU09cePH69jx47pww8/1KxZs3Tw4EF169bN5ZZGd78bAHgwAwAeKC4uzkgyjzzyiFvb3XvvvcbPz8+cPXvWGGPMnDlzjCQze/ZsZ50FCxYYSeb9999Ps/21a9dclqxMmDDBSDJxcXHm2rVr5uzZs2bRokUmICDA9OnTx1nv4sWLpmTJkqZbt24u2ycnJ5v69eubu+++O02bI0eOdKn76aefGklm/vz5xhhjYmJijLe3txk6dKhLvfPnz5vg4GDTq1evDONOSkoyFy5cMAEBAeadd95xlv/73/82ksy6deuy7Lsx/zu+mS2VKlVy1l+5cqWRZKZMmeLSzsKFC40kM2vWLGdZ7dq1TevWrdPs8/777zd33XVXpnFl9xju27fPSDKDBg1yqffjjz8aSWb8+PHOsq5du7r0JdWQIUNM8eLFM40nO5KTk821a9fM3//+d1OqVCmTkpJijDHmyJEjxsvLyzz66KOZbt+6dWsjyXz77bcu5S+++KKRZH788UeX8ueee87YbDZz4MCBbPejTp06pnv37u52zaxbty7L8ySzX0kkmcGDB2e5n+zWS+Xn52eaNWuWpjx1LFKX5ORk57pOnTqZChUqmPj4eJdthgwZYvz8/MyZM2eMMf/r83333edSb9GiRUaSiYqKMsa4990AwPNx5QpAoREdHa1169apR48eKl68uCTp4YcfVmBgYJpbA9Nz7tw5+fj4uCwZzYJ3s+DgYPn4+KhEiRLq1auXGjVqpI8//ti5fvPmzTpz5owGDBjg8hf7lJQUde7cWVu3bk1zpePRRx91+dyrVy95e3tr3bp1kqRVq1YpKSlJ/fv3d2nTz89PrVu31vr1653bXrhwQWPHjlXVqlXl7e0tb29vFStWTBcvXtS+ffuy1cfMzJs3T1u3bk2ztGzZ0qXed999J0lpJql4+OGHFRAQkK2/0t99993avXu3Bg0apFWrVikhISHDulkdw9T/3hzP3XffrZo1a2Y7nnPnzqlPnz768ssvnVe7suO7775T+/bt5XA45OXlJR8fH/3tb3/T6dOndeLECUnXr64kJydn6zbVEiVK6N57702zj1q1aunuu+92KR84cKCMMc4xyU4/7r77bn3zzTd68cUXtX79el2+fDnbfZWuP+uY3nnSq1cvt9rJbT169HD5Hhg2bJik67e0fvvtt/rLX/6iokWLuvzc3Xfffbpy5Yp++OEHl7ZuvnW5Xr16kuS86piT7wYAnosJLQB4pNKlS6to0aLO2+uy46OPPpIxRj179nS5ReyBBx7Qp59+qv3796tGjRqqWLGiJKW5pSowMFBbt26VdH2K9vRu8cnI2rVr5XA4dObMGc2aNUuLFy/W0KFD9f7770v637NXPXv2zLCNM2fOKCAgwPk5ODjYZb23t7dKlSrlvAUrtc0mTZqk216RIv/7+1nfvn317bff6pVXXlGTJk0UFBQkm82m++67z+1fkNNTs2bNdG+tcjgcLjMunj59Wt7e3mkmLLDZbAoODs7W7WXjxo1TQECA5s+fr/fff19eXl5q1aqVJk+enCaGrI5h6n9DQkLS7Cc0NDTNOZKefv36KSkpSR988IEeeughpaSkqEmTJnr99dfVoUOHDLfbsmWLOnbsqDZt2uiDDz5QhQoV5Ovrq2XLlukf//iHc1xSnxGrUKFClrGk14/Tp0+nO/NjaGioc312+zFt2jRVqFBBCxcu1OTJk+Xn56dOnTrpzTffVLVq1bKMr0qVKumeJ+lNYJEXKlasmO4Yv/XWW3r55Zcluf58nT59WklJSZo+fbqmT5+ebps3J6WlSpVy+Zz6/Gjq+ObkuwGA5yK5AuCRvLy81K5dO33zzTf67bffsvzFMiUlxTnxQo8ePdKt89FHH2nKlClq1KiRSpQoof/85z/OZ4xS95n6i9/PP//sVrz169d3zhbYoUMHderUSbNmzdKTTz6pJk2aONdNnz49w5nMypUr5/I5Li5O5cuXd35OSkrS6dOnnb+spbb5xRdfqFKlShnGFh8fr6+++koTJkxwmcAjMTExw/eE5ZZSpUopKSlJJ0+edPmF2hijuLi4DBPFG3l7e2vUqFEaNWqUzp07p7Vr12r8+PHq1KmTjh8/7jI7XlbHMPW/sbGxac6xP/74w3mMs/L444/r8ccf18WLF7Vx40ZNmDBB999/v3799dcMx2bBggXy8fHRV199JT8/P2f5smXLXOqlHqfffvsty+cPb5xkI1WpUqUUGxubpvyPP/6QJJc+ZtWPgIAAvfrqq3r11Vf1559/Oq9idevWTfv37880Nk/UoUMHvffee9q2bZtL0nfHHXekW79EiRLy8vJSv379MrySWLlyZbdiyMl3AwDPxW2BADzWuHHjZIzRU089patXr6ZZf+3aNf3nP/+RdP0Wud9++02DBw/WunXr0iy1a9fWvHnzlJSUJF9fX73wwgv6+eefNXnyZMvjttlseu+99+Tl5eX863eLFi1UvHhx7d27V40bN0538fX1dWnn008/dfm8aNEiJSUlqU2bNpKkTp06ydvbW4cPH86wzdR4jDFpZlz88MMP07wn6Oa/qlst9V1B8+fPdylfvHixLl686PIuIbvdnmUcxYsXV8+ePTV48GCdOXMmzcuPszqGqbfQ3RzP1q1btW/fPrfjCQgIUJcuXfTSSy/p6tWr+uWXXzKsa7PZ5O3tLS8vL2fZ5cuX9cknn7jU69ixo7y8vDRz5sxM952Rdu3aae/evdqxY4dL+bx582Sz2dJ9VUB2+lGuXDkNHDhQffr00YEDB3Tp0qUcxZefRo4cqaJFi2rw4ME6f/58lvWLFi2qtm3baufOnapXr166P3M3X6nKSk6+GwB4Lq5cAfBYzZs318yZMzVo0CA1atRIzz33nGrXrq1r165p586dmjVrlurUqaNu3bpp9uzZ8vb21vjx4523O93omWee0bBhw7RixQo9+OCDGjt2rPbv368XX3xRGzduVO/evRUeHq7ExEQdOXJEH374oby8vHL8jqBq1arp6aef1owZM7Rp0ya1bNlS06dP14ABA3TmzBn17NlTZcuW1cmTJ7V7926dPHkyzS/PS5Yskbe3tzp06KBffvlFr7zyiurXr+98PiU8PFx///vf9dJLL+nIkSPq3LmzSpQooT///FNbtmxxXmUICgpSq1at9Oabb6p06dIKDw/Xhg0bNHv2bOezaanq1KkjSZo1a5YCAwPl5+enypUru/0LY0ZSr+qNHTtWCQkJatGihfbs2aMJEyaoQYMG6tevn7Nu3bp1tWDBAi1cuFBVqlSRn5+f6tatq27duqlOnTpq3LixypQpo2PHjmnq1KmqVKlSmlvTsjqG1atX19NPP63p06erSJEi6tKli44ePapXXnlFYWFhGjlypEs8S5Ys0cyZM9WoUSMVKVJEjRs31lNPPSV/f3+1aNFCISEhiouLU2RkpBwOR6ZX4rp27aq3335bffv21dNPP63Tp0/rn//8Z5okODw8XOPHj9drr72my5cvq0+fPnI4HNq7d69OnTqV5e2rI0eO1Lx589S1a1f9/e9/V6VKlbRixQrNmDFDzz33nHNGy+z0o2nTprr//vtVr149lShRQvv27dMnn3yi5s2bF7j3aUnXr1B9/vnn6tOnj+rWravnnnvO+RLhEydOaPXq1ZLkMk37O++8o5YtW+qee+7Rc889p/DwcJ0/f16HDh3Sf/7zH+czbNlVrFgxt78bAHiwfJ1OAwCyYdeuXWbAgAGmYsWKxtfX1wQEBJgGDRqYv/3tb+bEiRPm5MmTxtfXN9NZzM6ePWv8/f3TzMi1fPly061bN1OuXDnj7e1tAgMDzV133WWef/55s3///ixjS52V7uTJk2nW/fnnn6ZYsWKmbdu2zrINGzaYrl27mpIlSxofHx9Tvnx507VrV/Pvf/87TZvbt2833bp1M8WKFTOBgYGmT58+5s8//0yzn2XLlpm2bduaoKAgY7fbTaVKlUzPnj3N2rVrnXV+++0389BDD5kSJUqYwMBA07lzZ/Pzzz+bSpUqmQEDBri0N3XqVFO5cmXj5eVlJJk5c+Zk2P/U2QK3bt2a7vr0Zti7fPmyGTt2rKlUqZLx8fExISEh5rnnnnPO8Jjq6NGjpmPHjiYwMNBl1sG33nrLREREmNKlSxtfX19TsWJF8+STT5qjR4/m6BgmJyebyZMnmzvvvNP4+PiY0qVLm8cee8wcP37cpd6ZM2dMz549TfHixY3NZnPObvfxxx+btm3bmnLlyhlfX18TGhpqevXqZfbs2ZPhcUv10UcfmerVqxu73W6qVKliIiMjzezZs40kEx0d7VJ33rx5pkmTJsbPz88UK1bMNGjQwGVsWrdubWrXrp3ufo4dO2b69u1rSpUqZXx8fEz16tXNm2++6TILXnb68eKLL5rGjRubEiVKOGMeOXKkOXXqVKb9TJ0578bz/EaDBw/Ol9kCUx0+fNgMHTrUVK9e3fj7+zt/jh5++GGzdOlS58yNqaKjo80TTzxhypcvb3x8fEyZMmVMRESEef311511MupzdHR0uj9X2fluAOD5bMbk08slAADpmjhxol599VWdPHky28/8wBXHEACQH3jmCgAAAAAsQHIFAAAAABbgtkAAAAAAsABXrgAAAADAAiRXAAAAAGABkisAAAAAsAAvEU5HSkqK/vjjDwUGBspms+V3OAAAAADyiTFG58+fV2hoqIoUyfzaFMlVOv744w+FhYXldxgAAAAAPMTx48dVoUKFTOuQXKUjMDBQ0vUDGBQUlM/RAAAAAMgvCQkJCgsLc+YImSG5SkfqrYBBQUEkVwAAAACy9bgQE1oAAAAAgAVIrgAAAADAAiRXAAAAAGABkisAAAAAsADJFQAAAABYgOQKAAAAACzAVOweLiU5RTHfx+h87HkFhgSq4j0VVcSLnBgAAADwNCRXHmzfkn1aOXylEn5LcJYFVQhS53c6q2aPmvkYGQAAAICbcQnEQ+1bsk+Lei5ySawkKeH3BC3quUj7luzLp8gAAAAApIfkygOlJKdo5fCVkkln5f+VrRyxUinJKXkaFwAAAICMkVx5oJjvY9JcsXJhpITjCYr5PibvggIAAACQKZIrD3Q+9ryl9QAAAADkPpIrDxQYEmhpPQAAAAC5j+TKA1W8p6KCKgRJtgwq2KSgsCBVvKdinsYFAAAAIGMkVx6oiFcRdX6n8/UPNydY//e589TOvO8KAAAA8CD8du6havaoqV5f9FJQ+SCX8qAKQer1RS/ecwUAAAB4GF4i7MFq9qip6g9WV8z3MTofe16BIYGqeE9FrlgBAAAAHojkysMV8Sqi8Dbh+R0GAAAAgCxwCQQAAAAALEByBQAAAAAWILkCAAAAAAuQXAEAAACABfI1uYqMjFSTJk0UGBiosmXLqnv37jpw4IBLHWOMJk6cqNDQUPn7+6tNmzb65Zdfsmx78eLFqlWrlux2u2rVqqWlS5fmVjcAAAAAIH+Tqw0bNmjw4MH64YcftGbNGiUlJaljx466ePGis86UKVP09ttv691339XWrVsVHBysDh066Pz58xm2GxUVpd69e6tfv37avXu3+vXrp169eunHH3/Mi24BAAAAuA3ZjDEmv4NIdfLkSZUtW1YbNmxQq1atZIxRaGioRowYobFjx0qSEhMTVa5cOU2ePFnPPPNMuu307t1bCQkJ+uabb5xlnTt3VokSJfT555+nqZ+YmKjExETn54SEBIWFhSk+Pl5BQUFp6gMAAAC4PSQkJMjhcGQrN/CoZ67i4+MlSSVLlpQkRUdHKy4uTh07dnTWsdvtat26tTZv3pxhO1FRUS7bSFKnTp0y3CYyMlIOh8O5hIWF3WpXAAAAANxmPCa5MsZo1KhRatmyperUqSNJiouLkySVK1fOpW65cuWc69ITFxfn1jbjxo1TfHy8czl+/PitdAUAAADAbcg7vwNINWTIEO3Zs0ebNm1Ks85ms7l8NsakKbuVbex2u+x2u5sRAwAAAMD/eMSVq6FDh2r58uVat26dKlSo4CwPDg6WpDRXnE6cOJHmytSNgoOD3d4GAAAAAG5FviZXxhgNGTJES5Ys0XfffafKlSu7rK9cubKCg4O1Zs0aZ9nVq1e1YcMGRUREZNhu8+bNXbaRpNWrV2e6DQAAAADciny9LXDw4MH67LPP9OWXXyowMNB5tcnhcMjf3182m00jRozQpEmTVK1aNVWrVk2TJk1S0aJF1bdvX2c7/fv3V/ny5RUZGSlJGj58uFq1aqXJkyfrwQcf1Jdffqm1a9eme8shAAAAAFghX5OrmTNnSpLatGnjUj5nzhwNHDhQkjRmzBhdvnxZgwYN0tmzZ9W0aVOtXr1agYGBzvoxMTEqUuR/F+EiIiK0YMECvfzyy3rllVd0xx13aOHChWratGmu9wkAAADA7cmj3nPlKdyZyx4AAABA4VVg33MFAAAAAAUVyRUAAAAAWIDkCgAAAAAsQHIFAAAAABYguQIAAAAAC5BcAQAAAIAFSK4AAAAAwAIkVwAAAABgAZIrAAAAALAAyRUAAAAAWIDkCgAAAAAsQHIFAAAAABYguQIAAAAAC5BcAQAAAIAFSK4AAAAAwAIkVwAAAABgAZIrAAAAALAAyRUAAAAAWIDkCgAAAAAsQHIFAAAAABYguQIAAAAAC5BcAQAAAIAFSK4AAAAAwAIkVwAAAABgAZIrAAAAALBAviZXGzduVLdu3RQaGiqbzaZly5a5rLfZbOkub775ZoZtzp07N91trly5ksu9AQAAAHA7y9fk6uLFi6pfv77efffddNfHxsa6LB999JFsNpseeuihTNsNCgpKs62fn19udAEAAAAAJEne+bnzLl26qEuXLhmuDw4Odvn85Zdfqm3btqpSpUqm7dpstjTbAgAAAEBuKjDPXP35559asWKFnnzyySzrXrhwQZUqVVKFChV0//33a+fOnZnWT0xMVEJCgssCAAAAAO4oMMnVxx9/rMDAQPXo0SPTejVq1NDcuXO1fPlyff755/Lz81OLFi108ODBDLeJjIyUw+FwLmFhYVaHDwAAAKCQsxljTH4HIV2/lW/p0qXq3r17uutr1KihDh06aPr06W61m5KSooYNG6pVq1aaNm1aunUSExOVmJjo/JyQkKCwsDDFx8crKCjIrf0BAAAAKDwSEhLkcDiylRvk6zNX2fX999/rwIEDWrhwodvbFilSRE2aNMn0ypXdbpfdbr+VEAEAAADc5grEbYGzZ89Wo0aNVL9+fbe3NcZo165dCgkJyYXIAAAAAOC6fL1ydeHCBR06dMj5OTo6Wrt27VLJkiVVsWJFSdcvw/373//WW2+9lW4b/fv3V/ny5RUZGSlJevXVV9WsWTNVq1ZNCQkJmjZtmnbt2qX33nsv9zsEAAAA4LaVr8nVtm3b1LZtW+fnUaNGSZIGDBiguXPnSpIWLFggY4z69OmTbhsxMTEqUuR/F+DOnTunp59+WnFxcXI4HGrQoIE2btyou+++O/c6AgAAAOC25zETWngSdx5aAwAAAFB4uZMbFIhnrgAAAADA05FcAQAAAIAFSK4AAAAAwAIkVwAAAABgAZIrAAAAALAAyRUAAAAAWIDkCgAAAAAsQHIFAAAAABYguQIAAAAAC5BcAQAAAIAFSK4AAAAAwAIkVwAAAABgAZIrAAAAALAAyRUAAAAAWIDkCgAAAAAsQHIFAAAAABYguQIAAAAAC5BcAQAAAIAFSK4AAAAAwAIkVwAAAABgAZIrAAAAALAAyRUAAAAAWIDkCgAAAAAsQHIFAAAAABYguQIAAAAAC5BcAQAAAIAF8jW52rhxo7p166bQ0FDZbDYtW7bMZf3AgQNls9lclmbNmmXZ7uLFi1WrVi3Z7XbVqlVLS5cuzaUeAAAAAMB1+ZpcXbx4UfXr19e7776bYZ3OnTsrNjbWuXz99deZthkVFaXevXurX79+2r17t/r166devXrpxx9/tDp8AAAAAHCyGWNMfgchSTabTUuXLlX37t2dZQMHDtS5c+fSXNHKTO/evZWQkKBvvvnGWda5c2eVKFFCn3/+ebbaSEhIkMPhUHx8vIKCgrK9bwAAAACFizu5gcc/c7V+/XqVLVtWd955p5566imdOHEi0/pRUVHq2LGjS1mnTp20efPmDLdJTExUQkKCywIAAAAA7vDo5KpLly769NNP9d133+mtt97S1q1bde+99yoxMTHDbeLi4lSuXDmXsnLlyikuLi7DbSIjI+VwOJxLWFiYZX0AAAAAcHvwzu8AMtO7d2/nv+vUqaPGjRurUqVKWrFihXr06JHhdjabzeWzMSZN2Y3GjRunUaNGOT8nJCSQYAEAAABwi0cnVzcLCQlRpUqVdPDgwQzrBAcHp7lKdeLEiTRXs25kt9tlt9stixMAAADA7SdHyVVKSooOHTqkEydOKCUlxWVdq1atLAksPadPn9bx48cVEhKSYZ3mzZtrzZo1GjlypLNs9erVioiIyLW4AAAAAMDt5OqHH35Q3759dezYMd080aDNZlNycnK227pw4YIOHTrk/BwdHa1du3apZMmSKlmypCZOnKiHHnpIISEhOnr0qMaPH6/SpUvrL3/5i3Ob/v37q3z58oqMjJQkDR8+XK1atdLkyZP14IMP6ssvv9TatWu1adMmd7sKAAAAANnmdnL17LPPqnHjxlqxYoVCQkIyfZYpK9u2bVPbtm2dn1OfexowYIBmzpypn376SfPmzdO5c+cUEhKitm3bauHChQoMDHRuExMToyJF/jcvR0REhBYsWKCXX35Zr7zyiu644w4tXLhQTZs2zXGcAAAAAJAVt99zFRAQoN27d6tq1aq5FVO+4z1XAAAAAKRcfs9V06ZNXW7lAwAAAABk87bAPXv2OP89dOhQPf/884qLi1PdunXl4+PjUrdevXrWRggAAAAABUC2bgssUqSIbDZbmgksnI383zp3J7TwVNwWCAAAAEByLzfI1pWr6OhoSwIDAAAAgMIqW8lVpUqVnP/euHGjIiIi5O3tumlSUpI2b97sUhcAAAAAbhduT2jRtm1bnTlzJk15fHy8y7TqAAAAAHA7cTu5Sn226manT59WQECAJUEBAAAAQEGT7ZcI9+jRQ9L1ySsGDhwou93uXJecnKw9e/YoIiLC+ggBAAAAoADIdnLlcDgkXb9yFRgYKH9/f+c6X19fNWvWTE899ZT1EQIAAABAAZDt5GrOnDmSpPDwcI0ePZpbAAEAAADgBtl6z9XthvdcAQAAAJBy4T1XN2rQoEG6E1rYbDb5+fmpatWqGjhwIDMHAgAAALituD1bYOfOnXXkyBEFBASobdu2atOmjYoVK6bDhw+rSZMmio2NVfv27fXll1/mRrwAAAAA4JHcvnJ16tQpPf/883rllVdcyl9//XUdO3ZMq1ev1oQJE/Taa6/pwQcftCxQAAAAAPBkbj9z5XA4tH37dlWtWtWl/NChQ2rUqJHi4+O1f/9+NWnSROfPn7c02LzCM1cAAAAAJPdyA7dvC/Tz89PmzZvTlG/evFl+fn6SpJSUFJf3YAEAAABAYef2bYFDhw7Vs88+q+3bt6tJkyay2WzasmWLPvzwQ40fP16StGrVKjVo0MDyYAEAAADAU+VoKvZPP/1U7777rg4cOCBJql69uoYOHaq+fftKki5fvuycPbAg4rZAAAAAAJJ7uQHvuUoHyRUAAAAAKZffc5Xq6tWrOnHihFJSUlzKK1asmNMmAQAAAKDAcju5OnjwoJ544ok0k1oYY2Sz2ZScnGxZcAAAAABQULidXA0cOFDe3t766quvFBISIpvNlhtxAQAAAECB4nZytWvXLm3fvl01atTIjXgAAAAAoEBy+z1XtWrV0qlTp3IjFgAAAAAosNxOriZPnqwxY8Zo/fr1On36tBISElwWAAAAALgduT0Ve5Ei1/Oxm5+1KkwTWjAVOwAAAAApl6diX7duXY4Du9nGjRv15ptvavv27YqNjdXSpUvVvXt3SdK1a9f08ssv6+uvv9aRI0fkcDjUvn17vfHGGwoNDc2wzblz5+rxxx9PU3758uUC+1JjAAAAAJ7P7eSqdevWlu384sWLql+/vh5//HE99NBDLusuXbqkHTt26JVXXlH9+vV19uxZjRgxQg888IC2bduWabtBQUE6cOCASxmJFQAAAIDclKOXCH///ff617/+pSNHjujf//63ypcvr08++USVK1dWy5Yts91Oly5d1KVLl3TXORwOrVmzxqVs+vTpuvvuuxUTE5Ppy4ptNpuCg4OzHQcAAAAA3Cq3J7RYvHixOnXqJH9/f+3YsUOJiYmSpPPnz2vSpEmWB3ij+Ph42Ww2FS9ePNN6Fy5cUKVKlVShQgXdf//92rlzZ6b1ExMTmZgDAAAAwC1xO7l6/fXX9f777+uDDz6Qj4+PszwiIkI7duywNLgbXblyRS+++KL69u2b6YNkNWrU0Ny5c7V8+XJ9/vnn8vPzU4sWLXTw4MEMt4mMjJTD4XAuYWFhudEFAAAAAIWY27MFFi1aVHv37lV4eLgCAwO1e/duValSRUeOHFGtWrV05cqVnAVis7lMaHGja9eu6eGHH1ZMTIzWr1/v1gx+KSkpatiwoVq1aqVp06alWycxMdF5BU66PiNIWFgYswUCAAAAt7lcnS0wJCREhw4dUnh4uEv5pk2bVKVKFXeby9K1a9fUq1cvRUdH67vvvnM72SlSpIiaNGmS6ZUru90uu91+q6ECAAAAuI25fVvgM888o+HDh+vHH3+UzWbTH3/8oU8//VSjR4/WoEGDLA0uNbE6ePCg1q5dq1KlSrndhjFGu3btUkhIiKWxAQAAAMCN3L5yNWbMGMXHx6tt27a6cuWKWrVqJbvdrtGjR2vIkCFutXXhwgUdOnTI+Tk6Olq7du1SyZIlFRoaqp49e2rHjh366quvlJycrLi4OElSyZIl5evrK0nq37+/ypcvr8jISEnSq6++qmbNmqlatWpKSEjQtGnTtGvXLr333nvudhUAAAAAss3tZ65SXbp0SXv37lVKSopq1aolu92u2NjYTKdIv9n69evVtm3bNOUDBgzQxIkTVbly5XS3W7dundq0aSNJatOmjcLDwzV37lxJ0siRI7VkyRLFxcXJ4XCoQYMGmjhxopo3b57tuNy5rxIAAABA4eVObpDj5Opmu3fvVsOGDZWcnGxFc/mK5AoAAACA5F5u4PYzVwAAAACAtEiuAAAAAMACJFcAAAAAYIFszxa4Z8+eTNcfOHDgloMBAAAAgIIq28nVXXfdJZvNpvTmv0gtt9lslgYHAAAAAAVFtpOr6Ojo3IwDAAAAAAq0bCdXlSpVys04AAAAAKBAY0ILAAAAALAAyRUAAAAAWIDkCgAAAAAsQHIFAAAAABZwO7m69957de7cuTTlCQkJuvfee62ICQAAAAAKHLeTq/Xr1+vq1atpyq9cuaLvv//ekqAAAAAAoKDJ9lTse/bscf577969iouLc35OTk7WypUrVb58eWujAwAAAIACItvJ1V133SWbzSabzZbu7X/+/v6aPn26pcEBAAAAQEGR7eQqOjpaxhhVqVJFW7ZsUZkyZZzrfH19VbZsWXl5eeVKkAAAAADg6bKdXFWqVEmSlJKSkmvBAAAAAEBBle3k6mZ79+5VTExMmsktHnjggVsOCgAAAAAKGreTqyNHjugvf/mLfvrpJ9lsNhljJEk2m03S9cktAAAAAOB24/ZU7MOHD1flypX1559/qmjRovrll1+0ceNGNW7cWOvXr8+FEAEAAADA87l95SoqKkrfffedypQpoyJFiqhIkSJq2bKlIiMjNWzYMO3cuTM34gQAAAAAj+b2lavk5GQVK1ZMklS6dGn98ccfkq5PeHHgwAFrowMAAACAAsLtK1d16tTRnj17VKVKFTVt2lRTpkyRr6+vZs2apSpVquRGjAAAAADg8dxOrl5++WVdvHhRkvT666/r/vvv1z333KNSpUpp4cKFlgcIAAAAAAWBzaRO93cLzpw5oxIlSjhnDCzoEhIS5HA4FB8fr6CgoPwOBwAAAEA+cSc3cPuZq1SHDh3SqlWrdPnyZZUsWTKnzQAAAABAoeB2cnX69Gm1a9dOd955p+677z7FxsZKkv7617/q+eefd6utjRs3qlu3bgoNDZXNZtOyZctc1htjNHHiRIWGhsrf319t2rTRL7/8kmW7ixcvVq1atWS321WrVi0tXbrUrbgAAAAAwF1uJ1cjR46Uj4+PYmJiVLRoUWd57969tXLlSrfaunjxourXr69333033fVTpkzR22+/rXfffVdbt25VcHCwOnTooPPnz2fYZlRUlHr37q1+/fpp9+7d6tevn3r16qUff/zRrdgAAAAAwB1uP3MVHBysVatWqX79+goMDNTu3btVpUoVRUdHq27durpw4ULOArHZtHTpUnXv3l3S9atWoaGhGjFihMaOHStJSkxMVLly5TR58mQ988wz6bbTu3dvJSQk6JtvvnGWde7cWSVKlNDnn3+erVh45goAAACAlMvPXF28eNHlilWqU6dOyW63u9tchqKjoxUXF6eOHTs6y+x2u1q3bq3NmzdnuF1UVJTLNpLUqVOnTLdJTExUQkKCywIAAAAA7nA7uWrVqpXmzZvn/Gyz2ZSSkqI333xTbdu2tSywuLg4SVK5cuVcysuVK+dcl9F27m4TGRkph8PhXMLCwm4hcgAAAAC3I7ffc/Xmm2+qTZs22rZtm65evaoxY8bol19+0ZkzZ/Tf//7X8gBvnt7dGJPllO/ubjNu3DiNGjXK+TkhIYEECwAAAIBb3L5yVatWLe3Zs0d33323OnTooIsXL6pHjx7auXOn7rjjDssCCw4OlqQ0V5xOnDiR5srUzdu5u43dbldQUJDLAgAAAADucPvKlXQ9gXn11VetjsVF5cqVFRwcrDVr1qhBgwaSpKtXr2rDhg2aPHlyhts1b95ca9as0ciRI51lq1evVkRERK7GCwAAAOD2lqPk6uzZs5o9e7b27dsnm82mmjVr6vHHH3f7ZcIXLlzQoUOHnJ+jo6O1a9culSxZUhUrVtSIESM0adIkVatWTdWqVdOkSZNUtGhR9e3b17lN//79Vb58eUVGRkqShg8frlatWmny5Ml68MEH9eWXX2rt2rXatGlTTroKAAAAANni9m2BGzZsUOXKlTVt2jSdPXtWZ86c0bRp01S5cmVt2LDBrba2bdumBg0aOK9MjRo1Sg0aNNDf/vY3SdKYMWM0YsQIDRo0SI0bN9bvv/+u1atXKzAw0NlGTEyM80XGkhQREaEFCxZozpw5qlevnubOnauFCxeqadOm7nYVAAAAALLN7fdc1alTRxEREZo5c6a8vLwkScnJyRo0aJD++9//6ueff86VQPMS77kCAAAAIOXye64OHz6s559/3plYSZKXl5dGjRqlw4cPux8tAAAAABQCbidXDRs21L59+9KU79u3T3fddZcVMQEAAABAgeP2hBbDhg3T8OHDdejQITVr1kyS9MMPP+i9997TG2+8oT179jjr1qtXz7pIAQAAAMCDuf3MVZEimV/sstlszpf2Jicn31Jw+YVnrgAAAABI7uUGbl+5io6OznFgAAAAAFBYuZ1cVapUKTfiAAAAAIACze0JLSTpk08+UYsWLRQaGqpjx45JkqZOnaovv/zS0uAAAAAAoKBwO7maOXOmRo0apfvuu0/nzp1zPldVvHhxTZ061er4AAAAAKBAcDu5mj59uj744AO99NJLLu+6aty4sX766SdLgwMAAACAgsLt5Co6OloNGjRIU26323Xx4kVLggIAAACAgsbt5Kpy5cratWtXmvJvvvlGtWrVsiImAAAAAChw3J4t8IUXXtDgwYN15coVGWO0ZcsWff7554qMjNSHH36YGzECAAAAgMdzO7l6/PHHlZSUpDFjxujSpUvq27evypcvr3feeUePPPJIbsQIAAAA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cdZhrPLP+fPnzc6dO83OnTuNJPP222+bnTt3OmeYy6uxSZ3eu127dmbHjh1m7dq1pkKFCkzvfYPMxuratWvmgQceMBUqVDC7du1y+X0jMTHR2QZjVbiRXOWD9957z1SqVMn4+vqahg0bOqcDR+6RlO4yZ84cZ52UlBQzYcIEExwcbOx2u2nVqpX56aefXNq5fPmyGTJkiClZsqTx9/c3999/v4mJiXGpc/r0afPoo4+awMBAExgYaB599FFz9uzZPOhl4XVzcsVYeY7//Oc/pk6dOsZut5saNWqYWbNmuaxnrDxDQkKCGT58uKlYsaLx8/MzVapUMS+99JLLL3yMVf5Zt25duv+PGjBggDEmb8fm2LFjpmvXrsbf39+ULFnSDBkyxFy5ciU3u1+gZDZW0dHRGf6+sW7dOmcbjFXhZjPGmLy7TgYAAAAAhRPPXAEAAACABUiuAAAAAMACJFcAAAAAYAGSKwAAAACwAMkVAAAAAFiA5AoAAAAALEByBQAAAAAWILkCAAAAAAuQXAEAbltt2rTRiBEj8juMdK1fv142m03nzp3LtF54eLimTp2aJzEBADJHcgUAyJG4uDgNHz5cVatWlZ+fn8qVK6eWLVvq/fff16VLl/I7vGxZsmSJXnvttRxv36ZNG9lsNtlsNtntdt15552aNGmSkpOTbzm2iIgIxcbGyuFwSJLmzp2r4sWLp6m3detWPf3007e8PwDArfPO7wAAAAXPkSNH1KJFCxUvXlyTJk1S3bp1lZSUpF9//VUfffSRQkND9cADD+R3mFkqWbLkLbfx1FNP6e9//7uuXLmir776SsOGDZOXl5fGjh17S+36+voqODg4y3plypS5pf0AAKzDlSsAgNsGDRokb29vbdu2Tb169VLNmjVVt25dPfTQQ1qxYoW6devmrBsfH6+nn35aZcuWVVBQkO69917t3r3buX7ixIm666679Mknnyg8PFwOh0OPPPKIzp8/76yTmJioYcOGqWzZsvLz81PLli21detW5/rUW+hWrVqlBg0ayN/fX/fee69OnDihb775RjVr1lRQUJD69OnjclXt5tsCExMTNWbMGIWFhclut6tatWqaPXt2pseiaNGiCg4OVnh4uIYMGaJ27dpp2bJlkqSzZ8+qf//+KlGihIoWLaouXbro4MGDzm2PHTumbt26qUSJEgoICFDt2rX1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",
      "text/plain": [
       "<Figure size 1000x400 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "positions_only = [p[0] for p in positions]\n",
    "\n",
    "plt.figure(figsize=(10,4))\n",
    "plt.scatter(positions_only, repeat_counts, color=\"purple\")\n",
    "plt.title(\"CAG Repeat Hotspots across HTT Gene\")\n",
    "plt.xlabel(\"Genomic Position\")\n",
    "plt.ylabel(\"Repeat Length\")\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 72,
   "id": "69257081-d165-4387-a8d1-2fb4917f7feb",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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du3UtWtX+fr6av/+/bpw4YJiY2Mdlp+enq769etbXu86derY/12xYkVJUmJioipXrixJioiIUPny5e01+VnfIUOG6IknntBHH32k1q1b629/+5uqVq1qr/X19XWYrlixohITEyVJycnJOn78uO666y6H5d91110Oh1Uud+DAAd15550OP7+mTZvm/01AzrhXB24mxhTZSxMocpGSkpLrPDc3N4fp7I1ITq78JsjCPiZ+4MABRUVF5TivQYMGio+P18qVK7Vu3To9+OCDat26tb744gtlZWWpYsWK2rRpk9Pzsi+PHDZsmFavXq3JkyerWrVq8vHxUffu3ZWeni5JCggI0A8//KBNmzZpzZo1euWVVzRmzBht377dfvXJ8uXLValSJYflF8RJo5cfisjeIF9+xYufn59DfX7Wd8yYMerZs6eWL1+ulStXavTo0Vq4cKG6du3q9JrZr2uu+DBfGe7yCnxXPhcAShICRS6u3AAVRa2rNmzYoJ9++kmDBw/OtSYwMFA9evRQjx491L17d7Vr106nT59WgwYNlJCQIHd3d/t5FVfaunWrevfubd+gpqSkOAUkd3d3tW7dWq1bt9bo0aNVunRpbdiwQbGxsfLy8tLRo0cVExNTUKt8zfKzvpJUo0YN1ahRQ4MHD9bDDz+sOXPm2Nc/L4GBgQoLC9M333zjcB7Jtm3b1KRJkxyfEx0d7XTPkO+++y5f6wMARY1AUUKlpaUpISFBmZmZOnnypFatWqXx48erY8eOevTRR3N8zttvv62KFSuqXr16KlWqlD7//HOFhoaqdOnSat26tZo2baouXbpowoQJuvXWW3X8+HGtWLFCXbp0UaNGjVStWjUtXrxYnTp1ks1m08svv+ywF2DZsmU6fPiw7r33XgUHB2vFihXKysrSrbfeqoCAAD3//PMaPHiwsrKydPfddys5OVnbtm2Tv7+/4uLirtdbJ0lXXd9atWpp2LBh6t69u6KiovT7779r+/bteuCBB/L9GsOGDdPo0aNVtWpV1atXT3PmzNHu3bv18ccf51j/zDPP6M0339SQIUP09NNPa+fOnQ5X7ABAcUagKKFWrVqlihUryt3dXcHBwapbt66mTp2quLg4p8Ms2fz9/TVhwgT98ssvcnNzU+PGjbVixQp7/YoVKzRq1Cj16dNHf/zxh0JDQ3XvvfcqJCRE0qVA0qdP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m7Rq1SrZbDYNGjRI48aNU/Xq1VW9enWNGzdOvr6+6tmzZ1G2DQAArlCkgeLkyZN65JFHdOLECQUFBalOnTpatWqVYmNjJUkvvPCCUlNT1a9fP505c0Z33HGH1qxZo4CAgKJsGwAAXMHl+1DMnz9fPXr0kJeXl8N4enq6Fi5cqEcffbRAG7SK+1AANznuQ4GbSRHeh8LlQOHm5qYTJ06oQoUKDuOnTp1ShQoVlJmZ6XrHhYhAAdzkCBS4mZSkG1vldtnm77//rqCgIFcXBwAAbgD5Poeifv36stlsstlsatWqldzd//fUzMxMxcfHq127doXSJAAAKN7yHSi6dOkiSdq9e7fatm0rf39/+zxPT09FRkbqgQceKPAGAQBA8ZfvQDF69GhJUmRkpHr06CFvb+9CawoAAJQsLl82GhcXJ+nSVR2JiYnKyspymF+5cuWC6QwAAJQYLgeKX375RX369HH6xs/skzWL21UeAACg8LkcKHr37i13d3ctW7ZMFStWzPOLugAAwM3B5UCxe/du7dy5U7fddlth9AMAAEogl+9DER0drT///LMwegEAACWUy4FiwoQJeuGFF7Rp0yadOnVKycnJDg8AAHDzcfnW26VKXcogV547UVxPyuTW28BNjvO8cDMpwltvu3wOxcaNG6+5MQAAcGNyOVDExMQURh8AAKAEczlQbNmyJc/599577zU3AwAASiaXA0Xz5s2dxi4/n6K4nUMBAAAKn8tXeZw5c8bhkZiYqFWrVqlx48Zas2ZNYfQIAACKOZf3UAQFBTmNxcbGysvLS4MHD9bOnTsLpDEAAFByuLyHIjfly5fXwYMHC2pxAACgBHF5D8WePXscpo0xOnHihN544w3VrVu3wBoDAAAlh8uBol69erLZbLryflh33nmnZs+eXWCNAQCAksPlQBEfH+8wXapUKZUvX17e3t4F1hQAAChZXA4UERERhdEHAAAowa7ppMzNmzerU6dOqlatmqpXr677779fW7duLejeAABACeFyoFiwYIFat24tX19fDRw4UAMGDJCPj49atWqlTz75pDB6BAAAxZzL3zZas2ZNPfXUUxo8eLDD+FtvvaX3339fBw4cKNAGreLbRoGbHN82iptJEX7bqMt7KA4fPqxOnTo5jd9///1OJ2wCAICbg8uBIjw8XOvXr3caX79+vcLDwwukKQAAULK4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QAEAACwjUAAAAMsIFAAAwDICBQAAsIxAAQAALCNQAAAAywgUAADAMgIFAACwjEABAAAsI1AAAADLCBQAAMAyAgUAALCMQAEAACwjUAAAAMsIFAAAwDICBQAAsIxAAQAALCNQAAAAywgUAADAMgIFAACwjEABAAAsI1AAAADLCBQAAMAyAgUAALCMQAEAACwjUAAAAMsIFAAAwDICBQAAsIxAAQAALCsRgeLdd99VVFSUvL291bBhQ23durWoWwIAAJcp9oHis88+06BBgzRq1Cjt2rVL99xzj+677z4dPXq0qFsDAAD/p9gHirfeekuPP/64nnjiCdWsWVNTpkxReHi4Zs6cWdStAQCA/+Ne1A3kJT09XTt37tSLL77oMN6mTRtt27Ytx+ekpaUpLS3NPp2UlCRJSk5OLtjmLhTs4oDirMA/PwAKRwF/VrM/+8aYq9YW60Dx559/KjMzUyEhIQ7jISEhSkhIyPE548eP19ixY53Gw8PDC6VH4GYQ9EZQUbcAID+CCuezeu7cOQVdZdnFOlBks9lsDtPGGKexbCNGjNCQIUPs01lZWTp9+rTKli2b63NQMiQnJys8PFzHjh1TYGBgUbcDIBd8Vm8cxhidO3dOYWFhV60t1oGiXLlycnNzc9obkZiY6LTXIpuXl5e8vLwcxkqXLl1YLaIIBAYG8kcKKAH4rN4YrrZnIluxPinT09NTDRs21Nq1ax3G165dq2bNmhVRVwAA4ErFeg+FJA0ZMkSPPPKIGjVqpKZNm2rWrFk6evSonnnmmaJuDQAA/J9iHyh69OihU6dO6dVXX9WJEydUu3ZtrVixQhEREUXdGq4zLy8vjR492umQFoDihc/qzclm8nMtCAAAQB6K9TkUAACgZCBQAAAAywgUAADAMgIFbgibNm2SzWbT2bNni7oVALgpESjgoHfv3rLZbHrjjTccxpcuXVri7zQaGRmpKVOmOI2PGTNG9erVs9fYbLZcH82bN89zvs1m05EjR67reuHmlf15zeky+n79+slms6l3794uLdNms2np0qUF0+Bljhw5IpvNpt27d19zXfPmzTVo0CB7jZXPamRkZIGv482u2F82iuvP29tbEyZM0NNPP63g4OACW256ero8PT0LbHmFYfv27crMzJQkbdu2TQ888IAOHjxov9vflevQuHFjPfXUU3ryySftY+XLl7++TeOmFh4eroULF+rtt9+Wj4+PJOnChQv69NNPVbly5SLurnCEh4frxIkT9unJkydr1apVWrdunX3s8s/qsWPH1KRJE61bt061atWSJLm5uV3fpm8C7KGAk9atWys0NFTjx4/Ps+7LL79UrVq15OXlpcjISL355psO8yMjI/X666+rd+/eCgoK0pNPPqm5c+eqdOnSWrZsmW699Vb5+vqqe/fuOn/+vObNm6fIyEgFBwfr2WeftW/YJWnBggVq1KiRAgICFBoaqp49eyoxMbHA1718+fIKDQ1VaGioypQpI0mqUKGCfaxy5cr2f4eGhsrNzc3e0+VjwPXSoEEDVa5cWYsXL7aPLV68WOHh4apfv75DbU576erVq6cxY8bY50tS165dHf4Xf+jQIXXu3FkhISHy9/dX48aNHTbe2c8dN26c+vTpo4CAAFWuXFmzZs2yz4+KipIk1a9f374H4Vq5ubk5fOb8/f3l7u7uMHb5ZzU75JctW9ZpDAWHQAEnbm5uGjdunKZNm6bff/89x5qdO3fqwQcf1EMPPaSffvpJY8aM0csvv6y5c+c61E2aNEm1a9fWzp079fLLL0uS/vrrL02dOlULFy7UqlWrtGnTJnXr1k0rVqzQihUr9NFHH2nWrFn64osv7MtJT0/Xa6+9ph9//FFLly5VfHy8y7tygRvVY489pjlz5tinZ8+erT59+ri8nO3bt0uS5syZoxMnTtinU1JS1L59e61bt067du1S27Zt1alTJx09etTh+W+++aYaNWqkXbt2qV+/furbt69+/vlnSdL3338vSVq3bp1OnDjhEIBwY+CQB3LUtWtX1atXT6NHj9aHH37oNP+tt95Sq1at7CGhRo0a2r9/vyZNmuSwoW/ZsqWef/55+/Q333yjixcvaubMmapataokqXv37vroo4908uRJ+fv7Kzo6Wi1atNDGjRvVo0cPSXL441ilShVNnTpVTZo0UUpKivz9/fO9XsOHD9dLL73kMJaenq7o6Oh8LwMobh555BGNGDHCfm7Bv//9by1cuFCbNm1yaTnZ/2svXbq0QkND7eN169ZV3bp17dOvv/66lixZoq+++koDBgywj7dv3179+vWTdOmz9vbbb2vTpk267bbbnPYSXE2zZs1UqpTj/3lTU1Pt5zuh+CFQIFcTJkxQy5YtNXToUKd5Bw4cUOfOnR3G7rrrLk2ZMkWZmZn23f6NGjVyeq6vr689TEhSSEiIIiMjHYJBSEiIwyGNXbt2acyYMdq9e7dOnz6trKwsSdLRo0ddCgPDhg1z2rMxdepUbdmyJd/LAIqbcuXKqUOHDpo3b56MMerQoYPKlStXYMs/f/68xo4dq2XLlun48ePKyMhQamqq0x6KOnXq2P9ts9kUGhp6zYcmP/vsM9WsWdNhrFevXte0LFwfBArk6t5771Xbtm01cuRIp42wMcbpqo+c7uLu5+fnNObh4eEwbbPZchzLDg3nz59XmzZt1KZNGy1YsEDly5fX0aNH1bZtW6Wnp7u0TuXKlVO1atUcxrLPlQBKsj59+tj3FsyYMSPHmlKlSjl9Ti9evHjVZQ8bNkyrV6/W5MmTVa1aNfn4+Kh79+5On7+8PseuCg8Pd/qsZp90iuKJQIE8vfHGG6pXr55q1KjhMB4dHa1vvvnGYWzbtm2qUaNGgZ+U+PPPP+vPP//UG2+8ofDwcEnSjh07CvQ1gJKuXbt29g1827Ztc6wpX768w9URycnJio+Pd6jx8PBwOCFakrZu3arevXura9euki6dU+Hq5dHZV1xcuWzcODgpE3m6/fbb1atXL02bNs1hfOjQoVq/fr1ee+01/fe//9W8efM0ffp0h/MlCkrlypXl6empadOm6fDhw/rqq6/02muvFfjrACWZm5ubDhw4oAMHDuQa6lu2bKmPPvpIW7du1d69exUXF+dUGxkZqfXr1yshIUFnzpyRJFWrVk2LFy/W7t279eOPP6pnz54u73moUKGCfHx8tGrVKp08eVJJSUnXtqIotggUuKrXXnvNaTdpgwYNtGjRIi1cuFC1a9fWK6+8oldffbVQrrwoX7685s6dq88//1zR0dF64403NHny5AJ/HaCkCwwMtN8zJScjRozQvffeq44dO6p9+/bq0qWLw/lM0qUrNdauXetw2enbb7+t4OBgNWvWTJ06dVLbtm3VoEEDl3pzd3fX1KlT9c9//lNhYWFO52Ch5OPrywEAgGXsoQAAAJYRKAAAgGUECgAAYBmBAgAAWEagAAAAlhEoAACAZQQKAABgGYECAABYRqAAAACWESgAAIBlBAoAAGAZgQIAAFj2/wGgB4ogR8GqNgAAAABJRU5ErkJggg==",
      "text/plain": [
       "<Figure size 600x500 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.figure(figsize=(6,5))\n",
    "plt.bar([\"Normal HTT\",\"Mutant HTT\"], [normal_max, mutant_max],\n",
    "        color=[\"green\",\"red\"])\n",
    "\n",
    "plt.axhline(y=36, linestyle=\"--\", color=\"black\", label=\"Disease Threshold\")\n",
    "\n",
    "plt.title(\"CAG Repeat Expansion Comparison\")\n",
    "plt.ylabel(\"Repeat Count\")\n",
    "plt.legend()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 74,
   "id": "5488e361-db14-4dfa-98ef-ade5855c60a4",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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4cIF8a9//Uv4+PgIDw8PMXr0aHHs2LEq6wpRccdmdHS0UCqVVe6y2rlzpxg7dqzw9fUVarVatGrVSowdO7ZK/D/88IOIi4sTGo1GREREiNdee61OPea7ublddRmj0Sieeuop0apVK+Hi4iJ69OghNmzYUOWuLCGE2LJli+jevbvQarUCgK28/7w7si5lrKkH8uq2mZCQIAYOHCh0Op0AYPe5qElBQYF46aWXRI8ePew+a1OnThV//PGHEEKI2267TWg0mqv2wj558mShUqlEZmambd7Zs2fF7NmzRbt27YSrq6vQarUiMjJSTJo0Saxfv97uzs7q1PRZrjRp0iShUqlEUlKSEKLiDtU5c+aI6OhooVarha+vr+jZs6d47rnnRElJid02lyxZIhYtWiTCwsKERqMR3bt3F5s2bao2hvvuu0+0atVKqNVqERAQIAYMGCBefvllu+Vee+01ERUVJbRarejQoYP44IMPqv0c1vUY4fId0gCEQqEQnp6eokuXLuKhhx4Se/bsqfE9q/xfMhgM4uGHHxZxcXHC09NTuLq6inbt2okFCxaI0tJS23p5eXni9ttvF97e3kKSJFvcVzsGV7s78sCBA2L8+PHC3d1deHh4iClTpohLly7ZrV/b84wQNZ8rqvs/LCsrE/PmzRORkZFCrVaLkJAQ8cgjj4j8/Hy75SIjI8XYsWOrlOuf51RyTpIQtbgthYiInEZKSgqio6OxdOlSPPXUU3KHQ0T1xLsjiYiIiGTAJIyIiIhIBrwcSURERCQD1oQRERERyYBJGBEREZEMmIQRERERyaDZd9ZqtVqRnp4ODw+POg9ZQURERFRX4nKnxKGhoVAoaq7vavZJWHp6ut2wLERERESNIS0t7apjojb7JKxyvLe0tLTrGl6EiKherFYgLa1iOjwcuMqvYiJqHoqKihAeHn7NMWebfRJWeQnS09OTSRgRNb7SUiAurmK6pARwc5M3HiJqNNdqBsWfZEREREQyYBJGREREJAMmYUREREQyYBJGREREJAMmYUREREQyYBJGREREJINm30UFEZGsVCrg0Uf/niYiuoxnBCIiR9JqgZUr5Y6CiJwQL0cSERERyYA1YUREjiQEkJNTMe3vD1yjB20iajmYhBEROZJeDwQGVkxz2CIiugIvRxIRERHJgEkYERERkQyYhBERERHJgG3CGkhqaipy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SRgRkQOpyvSYNb4bZo3vxmGLiMgO24QRETmY2lgmdwhE5IRYE0ZEREQkAyZhRERERDJgEkZEREQkAyZhRERERDJgEkZEREQkA94dSUTkQEKhQFpcH9s0EVElJmFERA5k0brg2zf+J3cYROSE+LOMiIiISAZMwoiIiIhkwCSMiMiBVGV6zJjUDzMm9eOwRURkh23CiIgcTFeYL3cIROSEWBNGREREJAMmYUREREQyYBJGREREJAMmYUREREQyYBJGREREJAPeHUlE5EBCoUBm2862aSKiSkzCiIgcyKJ1wRcrvpM7DCJyQvxZRkRERCQDp07CFi9ejN69e8PDwwOBgYG47bbbcOrUKbnDIiIiIrpuTp2E7dy5EzNnzsTevXuxefNmmM1mjBw5EqWlpXKHRkRUKypDGe6bdhPum3YTVIYyucMhIifi1G3CNm7caPd89erVCAwMxIEDBzBo0CCZoiIiqgMh4HXpom2aiKiSUydh/1RYWAgA8PX1rXEZo9EIo9Foe15UVOTwuIiIiIjqyqkvR15JCIG5c+fihhtuQOfOnWtcbvHixfDy8rI9wsPDGzFKIiIiotppMknYrFmzcOTIEXzxxRdXXW7+/PkoLCy0PdLS0hopQiIiIqLaaxKXI2fPno0ffvgBu3btQlhY2FWX1Wq10Gq1jRQZERERUf04dRImhMDs2bOxfv167NixA9HR0XKHRERERNQgnDoJmzlzJj7//HN8//338PDwQGZmJgDAy8sLrq6uMkdHRFQLkoTcyNa2aSKiSk6dhL333nsAgCFDhtjNX716Ne65557GD4iIqI7MLq745IOf5Q6DiJyQUydhgn3qEBERUTPVZO6OJCIiImpOmIQRETmQylCGux8ci7sfHMthi4jIjlNfjiQiavKEgN/5JNs0EVEl1oQRERERyYBJGBEREZEMmIQRERERyYBJGBEREZEMmIQRERERyYB3RxIROZIkoTColW2aiKgSkzAiIgcyu7ji4/9tkzsMInJCvBxJREREJAMmYUREREQyYBJGRORASqMBU2b9C1Nm/QtKo0HucIjIibBNGBGRA0lWK4JPH7NNExFVYk0YERERkQyYhBERERHJgEkYERERkQyYhBERERHJgEkYERERkQx4dyQRkYPpvXzkDoGInBCTMCIiBzK76vDfb/bKHQYROSFejiQiIiKSAZMwIiIi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jQkTJkCn09W6/EQkHyZhRNQo2rRpA0mSkJiYiNtuu63G5T7++GMcP34cKtXfpyer1YqPPvoIDz30EFq3bg1Jkuwa8QNATEwMAMDV1bVOcfn7+6N169Zo3bo1vvvuO7Ru3Rr9+vXD8OHDYbVaoVQqceDAASiVSrv13N3d7Z5Xl1xWzrNarVi0aBEmTpxYZRkXFxecP38eN998Mx5++GG89NJL8PX1xe+//477778fJpOpxtg9PDxw8OBB7NixA7/99hteeOEFLFy4EPHx8byTkqgJYBJGRI3C19cXo0aNwsqVKzFnzpwqtTsFBQVIS0vD/v37sWPHDvj6+tq9NmjQIBw7dgydO3fGiBEjsGLFCsyePbvGWqL68PHxwezZs/HUU0/h0KFD6N69OywWC7KysnDjjTfWuJ7ZbMb+/fvRp08fAMCpU6dQUFCA9u3bAwB69OiBU6dOoXXr1tWuv3//fpjNZrz55ptQKCqa6n799dd2y2g0GlgslirrqlQqDB8+HMOHD8eCBQvg7e2Nbdu2VZvwEZFzYcN8Imo07777LiwWC/r06YPvvvsOZ86cQWJiIt5++230798fH330Efr06YNBgwahc+fOtscNN9xge71yO2azGb169cJXX32FxMREnDp1Cp9++ilOnjxZpdaqLmbOnIlTp07hu+++Q9u2bXHXXXfh7rvvxrp165CcnIz4+HgsWbIEv/zyi20dtVqN2bNnY9++fTh48CDuvfde9OvXz5aUvfDCC/jkk0+wcOFCHD9+HImJifjqq6/w73//GwAQGxsLs9mMd955B+fOncP//vc/vP/++3ZxRUVFoaSkBFu3bkVOTg70ej1++uknvP3220hISMD58+fxySefwGq1ol27dvUuPxE1IrkbpRFRy5Keni5mzpxpaxTfqlUrccstt4hNmzYJPz8/8frrr1e73ptvvin8/f2F0Wi0bWfWrFkiOjpaqNVq4e7uLvr06SOWLl0qSktLaxVLdV1UCCHEgw8+KDp16iQsFosoLy8XL7zwgoiKihJqtVoEBweLCRMmiCNHjggh/m4w/91334mYmBih0WjETTfdJFJSUuy2uXHjRjFgwADh6uoqPD09RZ8+fcSqVatsry9btkyEhIQIV1dXMWrUKPHJJ58IAHZddzz88MPCz8/P1kXF7t27xeDBg4WPj49wdXUVcXFx4quvvqpV2YlIfpIQTtg1NBFRE7FmzRo8/vjjKCgokDsUImpieDmSiIiISAZMwoioWfrss8/g7u5e7aNTp05yh0dEBF6OJKJmqbi4GJcuXar2NbVajcjIyEaOiIjIHpMwIiIiIhnwciQRERGRDJiEEREREcmASRgRERGRDJiEEREREcmASRgRERGRDJiEEREREcmASRgRERGRDP4/xq0aDS4C2x4AAAAASUVORK5CYII=",
      "text/plain": [
       "<Figure size 700x500 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "import seaborn as sns\n",
    "\n",
    "plt.figure(figsize=(7,5))\n",
    "sns.histplot(patients[\"CAG_Repeats\"], bins=20, kde=True)\n",
    "plt.axvline(36, color=\"red\", linestyle=\"--\")\n",
    "plt.title(\"Simulated Patient CAG Repeat Distribution\")\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 76,
   "id": "f6812598-9326-45cf-b66f-105bca74dcad",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAmUAAAHUCAYAAAB78V9qAAAAOXRFWHRTb2Z0d2FyZQBNYXRwbG90bGliIHZlcnNpb24zLjkuMiwgaHR0cHM6Ly9tYXRwbG90bGliLm9yZy8hTgPZAAAACXBIWXMAAA9hAAAPYQGoP6dpAABOL0lEQVR4nO3dd3hUVeL/8c+kTXqAAKGFGAgoShNQDFUsIKKCqGAFXFBQFAHRBQvtx4qiIu4KYkH54qJmFVHWggbpTamKwGqUkgAJoaaQRmbO7w82swwpZEJIrsn79TzzPDtnzrn3nHvGzYdbztiMMUYAAACoVF6V3QEAAAAQygAAACyBUAYAAGABhDIAAAALIJQBAABYAKEMAADAAghlAAAAFkAoAwAAsABCGQAAgAUQyoBKNH/+fNlsNtfL399f9erVU48ePTR9+nSlpqYWajN58mTZbLZK6G3FWrlypdux8fb2VkREhO666y7t3r27Qvpw7bXX6tprr3W937dvn2w2m+bPn+/Rdnbt2qXJkydr3759hT4bMmSILrnkkgvqZ1mdfXxtNpvCwsJ07bXX6quvvvJ4W+vXr9fkyZN18uTJQp/NmTPH42MGVEeEMsAC3n//fW3YsEHx8fGaPXu22rZtq5deekktWrTQsmXL3OoOGzZMGzZsqKSeVrwXXnhBGzZs0IoVK/TXv/5V8fHx6ty5sw4ePFjhfalfv742bNigPn36eNRu165dmjJlSpGh7Pnnn9fixYvLqYeeu/POO7VhwwatW7dOs2fPVkpKim699VaPg9n69es1ZcoUQhlwAXwquwMApJYtW6pDhw6u93fccYfGjBmjLl26qH///kpISFBERIQkqVGjRmrUqFFldbXCNWvWTNdcc40kqVu3bqpRo4aGDh2q+fPn69lnny2yTVZWlgIDA8u9L3a73dWX8tK0adNy3Z6nIiIiXGPq1KmTYmNjFRMTo1mzZnkcPiuSMUY5OTkKCAio7K4A5YYzZYBFNW7cWK+++qoyMjL01ltvucqLuny5fPlyXXvttQoPD1dAQIAaN26sO+64Q1lZWa46eXl5mjZtmi677DLZ7XbVqVNHDz74oI4cOeK2rbi4OPXs2VP169dXQECAWrRoofHjx+vUqVNu9fbs2aO7775bDRo0kN1uV0REhK6//npt37690PZiY2MVFBSk4OBg9erVS9u2bSvzcSkIEPv373c7Hlu3btWdd96pmjVruoKOMUZz5sxR27ZtFRAQoJo1a+rOO+/Unj173LZpjNGMGTMUFRUlf39/tWvXTt98802hfRd3+fI///mP7rnnHkVERMhut6tx48YaNGiQcnNzNX/+fN11112SpB49erguFRZso6jLlzk5OZowYYKio6Pl5+enhg0bauTIkYXOQl1yySW65ZZbtHTpUrVr104BAQG67LLL9N5775Xl0Eo6ExLr1KnjOr7x8fHq27evGjVqJH9/f8XExGj48OE6evSoq83kyZP11FNPSZKio6NdY1y5cqUuueQS7dy5U6tWrXKVnz3e9PR0jRs3zm2so0ePLvR9s9lseuyxxzR37ly1aNFCdrtd//d//+e6BWDFihV65JFHVLt2bYWHh6t///46dOhQmY8DUBk4UwZY2M033yxvb2+tXr262Dr79u1Tnz591LVrV7333nuqUaOGDh48qKVLlyovL0+BgYFyOp3q27ev1qxZo6efflqdOnXS/v37NWnSJF177bXavHmz64xDQkKCbr75Zo0ePVpBQUH6z3/+o5deekk//vijli9f7tY3h8OhGTNmqHHjxjp69KjWr1/vFhxeeOEFPffcc3rwwQf13HPPKS8vTy+//LK6du2qH3/8UZdffrnHx+T333+XJNWpU8etvH///rr77rs1YsQI1x/04cOHa/78+Ro1apReeuklHT9+XFOnTlWnTp30008/uc4+TpkyRVOmTNHQoUN15513KikpSQ899JAcDocuvfTSEvvz008/qUuXLqpdu7amTp2qZs2aKTk5WUuWLFFeXp769OmjF154Qc8884xmz56tdu3aSSr+DJkxRv369dP333+vCRMmqGvXrvr55581adIkbdiwQRs2bJDdbnfb/5NPPqnx48crIiJC7777roYOHaqYmBh169bN4+N74sQJHTt2TM2aNZMk/fHHH4qNjdWwYcMUFhamffv2aebMmerSpYt27NghX19fDRs2TMePH9c//vEPffbZZ6pfv74k6fLLL9fixYt15513KiwsTHPmzJEkV/+zsrLUvXt3HThwQM8884xat26tnTt3auLEidqxY4eWLVvm9g+Qzz//XGvWrNHEiRNVr1491a1bV5s2bZJ05rJ+nz599OGHHyopKUlPPfWU7r//frfvLGB5BkClef/9940ks2nTpmLrREREmBYtWrjeT5o0yZz9n+6nn35qJJnt27cXu42PPvrISDKLFi1yK9+0aZORZObMmVNkO6fTaU6fPm1WrVplJJmffvrJGGPM0aNHjSQza9asYveZmJhofHx8zOOPP+5WnpGRYerVq2cGDBhQbFtjjFmxYoWRZOLi4szp06dNVlaWWb16tYmJiTHe3t6uvhQcj4kTJ7q137Bhg5FkXn31VbfypKQkExAQYJ5++mljjDEnTpww/v7+5vbbb3ert27dOiPJdO/e3VW2d+9eI8m8//77rrLrrrvO1KhRw6SmphY7lk8++cRIMitWrCj02eDBg01UVJTr/dKlS40kM2PGDLd6cXFxRpJ5++23XWVRUVHG39/f7N+/31WWnZ1tatWqZYYPH15sfwpIMo8++qg5ffq0ycvLM7t37za9e/c2kszs2bML1S/4Puzfv99IMl988YXrs5dfftlIMnv37i3U7oorrnA7jgWmT59uvLy8Cn3/C77TX3/9tVtfw8LCzPHjx93qFvw39Oijj7qVz5gxw0gyycnJ5z0OgFVw+RKwOGNMiZ+3bdtWfn5+evjhh/V///d/hS7NSdKXX36pGjVq6NZbb1V+fr7r1bZtW9WrV08rV6501d2zZ4/uvfde1atXT97e3vL19VX37t0lyfXUY61atdS0aVO9/PLLmjlzprZt2yan0+m2z2+//Vb5+fkaNGiQ2z79/f3VvXt3t32WZODAgfL19VVgYKC6desmh8OhTz/9VK1bt3ard8cddxQas81m0/333++2/3r16qlNmzau/W/YsEE5OTm677773Np36tRJUVFRJfYtKytLq1at0oABAwqduSurgjM7Q4YMcSu/6667FBQUpO+//96tvG3btmrcuLHrvb+/v5o3b+66/Hg+c+bMka+vr/z8/NSiRQutX79eU6dO1aOPPipJSk1N1YgRIxQZGSkfHx/5+vq6jsuFPgX75ZdfqmXLlmrbtq3bHPXq1ct1+fNs1113nWrWrFnktm677Ta39wXfj9IeB8AKuHwJWNipU6d07NgxtWrVqtg6TZs21bJlyzRjxgyNHDlSp06dUpMmTTRq1Cg98cQTkqTDhw/r5MmT8vPzK3IbBfcHZWZmqmvXrvL399e0adPUvHlzBQYGKikpSf3791d2drakM/f3fP/995o6dapmzJihJ598UrVq1dJ9992n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U5WcnKyWLVvq66+/VlRUlCQpOTnZbc2yIUOGKCMjQ2+88YaefPJJ1ahRQ9ddd51eeumlyhoCAABAubCZanbdLz09XWFhYUpLS1NoaGhldwcAAFRxpc0elf70JQAAAAhlAAAAlkAoAwAAsABCGQAAgAUQygAAACyAUAYAAGABhDIAAAALIJQBAABYAKEMAADAAghlAAAAFkAoAwAAsABCGQAAgAUQygAAACyAUAYAAGABhDIAAAALIJQBAABYAKEMAADAAghlAAAAFkAoAwAAsABCGQAAgAUQygAAACyAUAYAAGABhDIAAAALIJQBAABYAKEMAADAAghlAAAAFkAoAwAAsABCGQAAgAUQygAAACyAUAYAAGABhDIAAAALIJQBAABYAKEMAADAAghlAAAAFkAoAwAAsABCGQAAgAUQygAAACyAUAYAAGABhDIAAAALIJQBAABYAKEMAADAAghlAAAAFkAoAwAAsABCGQAAgAUQygAAACyAUAYAAGABhDIAAAALIJQBAABYAKEMAADAAghlAAAAFkAoAwAAsABCGQAAgAX4lKbSzz//XOoNtm7dusydAQAAqK5KFcratm0rm80mY0yRnxd8ZrPZ5HA4yrWDAAAA1UGpQtnevXsvdj8AAACqtVKFsqioqIvdDwAAgGqtTDf6f/DBB+rcubMaNGig/fv3S5JmzZqlL774olw7BwAAUF14HMrefPNNjR07VjfffLNOnjzpuoesRo0amjVrVnn3DwAAoFrwOJT94x//0DvvvKNnn31W3t7ervIOHTpox44d5do5AACA6sLjULZ3715deeWVhcrtdrtOnTrlcQfmzJmj6Oho+fv7q3379lqzZk2J9XNzc/Xss88qKipKdrtdTZs21XvvvefxfgEAAKykVDf6ny06Olrbt28vdPP/N998o8svv9yjbcXFxWn06NGaM2eOOnfurLfeeku9e/fWrl271Lhx4yLbDBgwQIcPH9a8efMUExOj1NRU5efnezoMAAAAS/E4lD311FMaOXKkcnJyZIzRjz/+qI8++kjTp0/Xu+++69G2Zs6cqaFDh2rYsGGSzjws8O233+rNN9/U9OnTC9VfunSpVq1apT179qhWrVqSpEsuucTTIQAAAFiOx6HswQcfVH5+vp5++mllZWXp3nvvVcOGDfX666/r7rvvLvV28vLytGXLFo0fP96t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JHXv3l133XXXxegjAABAledxKNu9e7cGDx4sSfLx8VF2draCg4M1depUvfTSS+XeQQAAgOrA41AWFBSk3NxcSVKDBg30xx9/uD47evRo+fUMAACgGvH4nrJrrrlG69at0+WXX64+ffroySef1I4dO/TZZ5/pmmuuuRh9BAAAqPI8DmUzZ85UZmamJGny5MnKzMxUXFycYmJi9Nprr5V7BwEAAKoD1ikDAAC4iEqbPTy+p0ySTp48qXfffVcTJkzQ8ePHJUlbt27VwYMHy9ZbAACAas7jy5c///yzbrjhBoWFhWnfvn166KGHVKtWLS1evFj79+/XggULLkY/AQAAqjSPz5SNHTtWQ4YMUUJCgvz9/V3lvXv31urVq8u1cwAAANWFx6Fs06ZNGj58eKHyhg0bKiUlpVw6BQAAUN14HMr8/f2Vnp5eqPzXX39VnTp1yqVTAAAA1Y3Hoaxv376aOnWqTp8+LUmy2WxKTEzU+PHjdccdd5R7BwEAAKoDj0PZK6+8oiNHjqhu3brKzs5W9+7dFRMTo5CQEP3tb3+7GH0EAACo8jx++jI0NFRr167V8uXLtXXrVjmdTrVr10433HDDxegfAABAteBxKCtw3XXX6brrrpN0Zt0yAAAAlJ3Hly9feuklxcXFud4PGDBA4eHhatiwoX766ady7RwAAEB14XEoe+uttxQZGSlJio+PV3x8vL755hv17t1bTz31VLl3EAAAoDrwOJQlJye7QtmXX36pAQMGqGfPnnr66ae1adMmjzswZ84cRUdHy9/fX+3bt9eaNWtK1W7dunXy8fFR27ZtPd4nAACA1XgcymrWrKmkpCRJ0tKlS103+Btj5HA4PNpWXFycRo8erWeffVbbtm1T165d1bt3byUmJpbYLi0tTYMGDdL111/vafcBAAAsyeNQ1r9/f91777268cYbdezYMfXu3VuStH37dsXExHi0rZkzZ2ro0KEaNmyYWrRooVmzZikyMlJvvvlmie2GDx+ue++9V7GxsZ52HwAAwJI8DmWvvfaaHnvsMV1++eWKj49XcHCwpDOXNR999NFSbycvL09btmxRz5493cp79uyp9evXF9vu/fff1x9//KFJkyaVaj+5ublKT093ewEAAFiNx0ti+Pr6aty4cYXKR48e7dF2jh49KofDoYiICLfyiIiIYn9DMyEhQePHj9eaNWvk41O6rk+fPl1TpkzxqG8AAAAVrVTJZsmSJerdu7d8fX21ZMmSEuvedtttHnXAZrO5vTfGFCqTJIfDoXvvvVdTpkxR8+bNS739CRMmaOzYsa736enprgcVAAAArKJUoaxfv35KSUlR3bp11a9fv2Lr2Wy2Ut/sX7t2bXl7exc6K5aamlro7JkkZWRkaPPmzdq2bZsee+wxSZLT6ZQxRj4+Pvruu+9ci9mezW63y263l6pPAAAAlaVUoczpdBb5vy+En5+f2rdvr/j4eN1+++2u8vj4ePXt27dQ/dDQUO3YscOtbM6cOVq+fLk+/fRTRUdHl0u/AAAAKkOZf2apPIwdO1YPPPCAOnTooNjYWL399ttKTEzUiBEjJJ259Hjw4EEtWLBAXl5eatmypVv7unXryt/fv1A5AADAn41HoczpdGr+/Pn67LPPtG/fPtlsNkVHR+vOO+/UAw88UOS9YCUZOHCgjh07pqlTpyo5OVktW7bU119/raioKElnnug835plAAAAVYHNGGNKU9EYo1tvvVVff/212rRpo8suu0zGGO3evVs7duzQbbfdps8///wid/fCpaenKywsTGlpaQoNDa3s7gAAgCqutNmj1GfK5s+fr9WrV+v7779Xjx493D5bvny5+vXrpwULFmjQoEFl7zUAAEA1VerFYz/66CM988wzhQKZJF133XUaP368Fi5cWK6dAwAAqC5KHcp+/vln3XTTTcV+3rt3b/3000/l0ikAAIDqptSh7Pjx40WuH1YgIiJCJ06cKJdOAQAAVDelDmUOh6PEnzby9vZWfn5+uXQKAACguin1jf7GGA0ZMqTY1fFzc3PLrVMAAADVTalD2eDBg89bhycvAQAAyqbUoez999+/mP0AAACo1kp9TxkAAAAuHkIZAACABRDKAAAALIBQBgAAYAGEMgAAAAsglAEAAFgAoQwAAMACCGUAAAAWQCgDAACwAEIZAACABRDKAAAALIBQBgAAYAGEMgAAAAsglAEAAFgAoQwAAMACCGUAAAAWQCgDAACwAEIZAACABRDKAAAALIBQBgAAYAGEMgAAAAsglAEAAFgAoQwAAMACCGUAAAAWQCgDAACwAEIZAACABRDKAAAALIBQBgAAYAGEMgAAAAsglAEAAFgAoQwAAMACCGUAAAAWQCgDAACwAEIZAACABRDKAAAALIBQBgAAYAGEMgAAAAsglAEAAFgAoQwAAMACCGUAAAAWQCgDAACwAEIZAACABRDKAAAALIBQBgAAYAGEMgAAAAsglAEAAFgAoQwAAMACCGUAAAAWQCgDAACwAEIZAACABRDKAAAALKDSQ9mcOXMUHR0tf39/tW/fXmvWrCm27meffaYbb7xRderUUWhoqGJjY/Xtt99WYG8BAAAujkoNZXFxcRo9erSeffZZbdu2TV27dlXv3r2VmJhYZP3Vq1frxhtv1Ndff60tW7aoR48euvXWW7Vt27YK7jkAAED5shljTGXtvGPHjmrXrp3efPNNV1mLFi3Ur18/TZ8+vVTbuOKKKzRw4EBNnDixVPXT09MVFhamtLQ0hYaGlqnfAAAApVXa7FFpZ8ry8vK0ZcsW9ezZ0628Z8+eWr9+fam24XQ6lZGRoVq1ahVbJzc3V+np6W4vAAAAq6m0UHb06FE5HA5FRES4lUdERCglJaVU23j11Vd16tQpDRgwoNg606dPV1hYmOsVGRl5Qf0GAAC4GCr9Rn+bzeb23hhTqKwoH330kSZPnqy4uDjVrVu32HoTJkxQWlqa65WUlHTBfQYAAChvPpW149q1a8vb27vQWbHU1NRCZ8/OFRcXp6FDh+qTTz7RDTfcUGJdu90uu91+wf0FAAC4mCrtTJmfn5/at2+v+Ph4t/L4+Hh16tSp2HYfffSRhgwZog8//FB9+vS52N0EAACoEJV2pkySxo4dqwceeEAdOnRQbGys3n77bSUmJmrEiBGSzlx6PHjwoBYsWCDpTCAbNGiQXn/9dV1zzTWus2wBAQEKCwurtHEAAABcqEoNZQMHDtSxY8c0depUJScnq2XLlvr6668VFRUlSUpOTnZbs+ytt95Sfn6+Ro4cqZEjR7rKBw8erPnz51d09wEAAMpNpa5TVhlYpwwAAFQky69TBgAAgP8hlAEAAFgAoQwAAMACCGUAAAAWQCgDAACwAEIZAACABRDKAAAALIBQBgAAYAGEMgAAAAsglAEAAFgAoQwAAMACCGUAAAAWQCgDAACwAEIZAACABRDKAAAALIBQBgAAYAGEMgAAAAsglAEAAFgAoQwAAMACCGUAAAAWQCgDAACwAEIZAACABRDKAAAALIBQBgAAYAGEMgAAAAsglAEAAFgAoQwAAMACCGUAAAAWQCgDAACwAEIZAACABRDKAAAALIBQBgAAYAGEMgAAAAsglAEAAFgAoQwAAMACCGUAAAAWQCgDAACwAEIZAACABRDKAAAALIBQBgAAYAGEMgAAAAsglAEAAFgAoQwAAMACCGUAAAAWQCgDAACwAEIZAACABRDKAAAALIBQBgAAYAGEMgAAAAsglAEAAFgAoQwAAMACCGUAAAAWQCgDAACwAEIZAACABRDKAAAALIBQBgAAYAGEMgAAAAsglAEAAFiAT2V3oKrJPJWnSV/uUuKJLDWuGagpt1yu4CC/EtscT8/WsAVblJyRo/oh/np3UHvVCg0otn5W1mm9suw31z7G3dBcgYG+Je6jLG3y8hxavP2gDp7MVsMaAbq9bUP5+XlXer/SM3M1btHPSjqZrcgaAXrljtYKDbaXW31JcjqN9h07pYycfIX4++iS8CB5edlKbOPp3Ht6fMvSpizfx5MZORr50TYdTMtRwzB/zb7nStUI8S/XfuXk5OvddXuVdCJLkTUDNaxztPz9S/6/o7Icr/x8p9b9cVRHMnJVJ8Suzk1ry8en+H+Lelq/rGMBgKLYjDGmMjswZ84cvfzyy0pOTtYVV1yhWbNmqWvXrsXWX7VqlcaOHaudO3eqQYMGevrppzVixIhS7y89PV1hYWFKS0tTaGhoeQzB5f55G7U24Vih8i7NwvXPodcU2abbjOVKPJ5dqLxxrQCtfvq6QuUjP9yib3akyHnWrHnZpN6t6mn2ve2L3EdZ2sxekaB5a/YqPSdfTmPkZbMp1N9HQ7tGa2SPZpXWr36z12p7Ulqh8raRYfp8ZJcLri9JvxxM06KtB/R7aqZyTztl9/VSTN1g3dGukVo2DCuyjadz7+nxLUubsnwfe762Sr8dzixU3jwiWN+N6V4u/Zr4xQ7FbUpSbv7/Jt7uY9PAqyI1tW+rctmHJH2x/aDeWvWHDqXlKN9h5ONtU4Mwfw3v3lR92za84PplHQuA6qe02aNSL1/GxcVp9OjRevbZZ7Vt2zZ17dpVvXv3VmJiYpH19+7dq5tvvlldu3bVtm3b9Mwzz2jUqFFatGhRBfe8sOL+AErS2oRjun/exkLlxQUySUo8nq1uM5a7lY38cIu++tk9xEiS00hf/ZyikR9uKbSdsrSZvSJBry9L0Mns0/LxtinQ11s+3jadzD6t15claPaKhErpV3EBS5K2J6Wp3+y1F1RfOhPI/v59gnYcSFONAD9dUjtINQL8tOPAmfJfDhbenqdz7+nxLUubsnwfiwtkkvTb4Uz1fG3VBfdr4hc79MGGROXmG3nZJG/bmSCem2/0wYZETfxixwXvQzoTsKZ9uUuJx7IU6OutOsF+CvT1VuKxLE37cpe+2H7wguqXdSwAUJJKDWUzZ87U0KFDNWzYMLVo0UKzZs1SZGSk3nzzzSLrz507V40bN9asWbPUokULDRs2TH/5y1/0yiuvVHDP3WWeyiv2D2CBtQnHlHkqz/X+eHp2sYGsQOLxbB1PP1MnK+u0vtmR4vrMy/a/V4FvdqQoK+u0631Z2uTlOTRvzV7lO40C/bzl5+0lLy+b/Ly9FOjnrXyn0by1e5WX56jQfqVn5hYbsApsT0pTemZumepLZy5ZLtp6QMdP5SmmbrCC/X3k7WVTsL+PYuoG6/ipPH229aCcZyVJT+fe0+NbljZl+T6ezMgpNpAV+O1wpk5m5JS5Xzk5+YrblCSjMwHG28smLy+bvL1s8rZJRlLcpiTl5OSXeR/SmUuQb636Q9l5DtUNtSvAz0deXl4K8PNR3VC7svMcenvVHuXnO8tUv6xjAYDzqbRQlpeXpy1btqhnz55u5T179tT69euLbLNhw4ZC9Xv16qXNmzfr9OnTRbbJzc1Venq626u8Tfpyl8f1hi0ofCaoKAX1Xln2m+us0rm3NhW8d5oz9QqUpc3i7QeVnpMvPx8vedls57Sxyc/HS+nZ+Vr83zMHFdWvcYt+VmkU1PO0viTtO3ZKv6dmqn5YgGznjN1ms6l+WIASUjO079gpV7mnc+/p8S1Lm7J8H0d+tK1Ubc6u52m/3l2313VW6dz787y8bK6zTO+u21vmfUjSuj+O6lBajkIDfGWzuf9fnM3mpdAAXx1My9a6P46WqX5ZxwIA51Npoezo0aNyOByKiIhwK4+IiFBKSkqRbVJSUoqsn5+fr6NHjxbZZvr06QoLC3O9IiMjy2cAZ0k8keVxveSzzjiUpKBeWfZRljYHT2bLaYx8zvkDWMDbZpPTGB08mV2h/Uo6WfJZxXPreVpfkjJy8pV72qmAYm4eD/DzVu5ppzLOOvvh6Vg8Pb5laVOmeU8r3ffx7Hqe9ivpv/sr7nGJgvKkC/g+StKRjFzlO4zsxdyg7+fjpXyH0ZGM3DLVL+tYAOB8Kn1JjHPPSBhjCpWdr35R5QUmTJigtLQ01yspKekCe1xY45qBHterf56n2c6tV5Z9lKVNwxoB8rLZlF/M8x+O/95k3bBGQIX2K7JG8U+jnq2gnqf1JSnE30d2Xy9ln3Up7GzZeQ7Zfb0UctaTdZ6OxdPjW5Y2ZZr3sNJ9H8+u52m/Iv+7v+KeLCooj7yA76Mk1Qmxy8fbptyzLjeeLS/fKR9vm+qE2MtUv6xjAYDzqbRQVrt2bXl7exc6K5aamlrobFiBevXqFVnfx8dH4eHhRbax2+0KDQ11e5W3Kbdc7nG9dwcV/XThuQrqjbuhudulvbOdfSlw3A3NXeVlaXN724YK9fdRXr5TznP+EDqNUV6+U6EBPrr9v0+jVVS/XrmjtUqjoJ6n9SXpkvAgxdQNVnJats59KNkYo+S0bDWrG6JLwoNc5Z7OvafHtyxtyvJ9nH3PlaVqc3Y9T/s1rHO07D42OY3c7suTzrx3mjNPLg7rHF3mfUhS56a11SDMX+nZp2WMe9Ayxqn07NNqGBagzk1rl6l+WccCAOdTaaHMz89P7du3V3x8vFt5fHy8OnXqVGSb2NjYQvW/++47dejQQb6+Ja9tdTEFB/mpS7OiQ2GBLs3C3daHqhUaoMa1Sj6b07hWgGu9ssBAX/VuVc/1mdP871Wgd6t6bmt8laWNn5+3hnaNlo+XTVl5DuU5nHI4jfIcTmXlOeTjZdPQLtGu9aEqql+hwXa1jSx6OYoCbSPDXOuPeVpfOnMv0B3tGqlWkJ9+T81UZk6+HE6jzJx8/Z6aqVpBfurfrqHbPUSezr2nx7csbcryfawR4q/mEcEltmkeEey2Xpmn/fL399HAqyJlk+QwksNp5HQaOZxGDnPmkt/AqyLd1vgqy/Hy8fHS8O5NFeDnrdT0XGXn5cvhdCo7L1+p6bkK8PPWw92buNYf87R+WccCAOdTqeuUxcXF6YEHHtDcuXMVGxurt99+W++884527typqKgoTZgwQQcPHtSCBQsknVkSo2XLlho+fLgeeughbdiwQSNGjNBHH32kO+64o1T7ZJ2y0rUpcl2oAB8N7VI91ylrVjdE/ds1vLjrlJVwfMvSplLXKSuhX+W2Ttl5jldR6441DAvQw92blHqdspLql3UsAKqf0mYPSyweO2PGDCUnJ6tly5Z67bXX1K1bN0nSkCFDtG/fPq1cudJVf9WqVRozZoxr8di//vWvllk8VmJFf1b0Z0V/VvRnRX8A7v40oayiXexQBgAAcLY/xYr+AAAAOINQBgAAYAGEMgAAAAsglAEAAFgAoQwAAMACCGUAAAAWQCgDAACwAEIZAACABRDKAAAALIBQBgAAYAHV7gfaCn5VKj09vZJ7AgAAqoOCzHG+X7asdqEsIyNDkhQZGVnJPQEAANVJRkaGwsLCiv282v0gudPp1KFDhxQSEiKbzXZR9pGenq7IyEglJSVVux89r85jl6r3+Bk7Y2fs1Ut1Hr+nYzfGKCMjQw0aNJCXV/F3jlW7M2VeXl5q1KhRhewrNDS02n1RC1TnsUvVe/yMnbFXN9V57FL1Hr8nYy/pDFkBbvQHAACwAEIZAACABRDKLgK73a5JkybJbrdXdlcqXHUeu1S9x8/YGXt1U53HLlXv8V+ssVe7G/0BAACsiDNlAAAAFkAoAwAAsABCGQAAgAUQygAAACyAUHYBVq9erVtvvVUNGjSQzWbT559/7vb5kCFDZLPZ3F7XXHNN5XS2HE2fPl1XXXWVQkJCVLduXfXr10+//vqrWx1jjCZPnqwGDRooICBA1157rXbu3FlJPS5fpRl/VZ37N998U61bt3YtmBgbG6tvvvnG9XlVnvfzjb2qznlRpk+fLpvNptGjR7vKqvLcn62osVfluZ88eXKhsdWrV8/1eVWe9/ON/WLMO6HsApw6dUpt2rTRG2+8UWydm266ScnJya7X119/XYE9vDhWrVqlkSNHauPGjYqPj1d+fr569uypU6dOuerMmDFDM2fO1BtvvKFNmzapXr16uvHGG12/PfpnVprxS1Vz7hs1aqQXX3xRmzdv1ubNm3Xdddepb9++rv8Trsrzfr6xS1Vzzs+1adMmvf3222rdurVbeVWe+wLFjV2q2nN/xRVXuI1tx44drs+q+ryXNHbpIsy7QbmQZBYvXuxWNnjwYNO3b99K6U9FSk1NNZLMqlWrjDHGOJ1OU69ePfPiiy+66uTk5JiwsDAzd+7cyurmRXPu+I2pPnNvjDE1a9Y07777brWbd2P+N3ZjqsecZ2RkmGbNmpn4+HjTvXt388QTTxhjqsd/88WN3ZiqPfeTJk0ybdq0KfKzqj7vJY3dmIsz75wpu8hWrlypunXrqnnz5nrooYeUmppa2V0qd2lpaZKkWrVqSZL27t2rlJQU9ezZ01XHbrere/fuWr9+faX08WI6d/wFqvrcOxwOffzxxzp16pRiY2Or1byfO/YCVX3OR44cqT59+uiGG25wK68Oc1/c2AtU5blPSEhQgwYNFB0drbvvvlt79uyRVD3mvbixFyjvea92P0hekXr37q277rpLUVFR2rt3r55//nldd9112rJlS5VZAdkYo7Fjx6pLly5q2bKlJCklJUWSFBER4VY3IiJC+/fvr/A+XkxFjV+q2nO/Y8cOxcbGKicnR8HBwVq8eLEuv/xy1/8JV+V5L27sUtWec0n6+OOPtXXrVm3atKnQZ1X9v/mSxi5V7bnv2LGjFixYoObNm+vw4cOaNm2aOnXqpJ07d1b5eS9p7OHh4Rdn3sv1vFs1piIuX57r0KFDxtfX1yxatKhiOlUBHn30URMVFWWSkpJcZevWrTOSzKFDh9zqDhs2zPTq1auiu3hRFTX+olSluc/NzTUJCQlm06ZNZvz48aZ27dpm586d1WLeixt7UarSnCcmJpq6deua7du3u8rOvoRXlef+fGMvSlWa+3NlZmaaiIgI8+qrr1bpeS/K2WMvSnnMO5cvK1D9+vUVFRWlhISEyu5KuXj88ce1ZMkSrVixQo0aNXKVFzydUvCvqAKpqamF/kX1Z1bc+ItSlebez89PMTEx6tChg6ZPn642bdro9ddfrxbzXtzYi1KV5nzLli1KTU1V+/bt5ePjIx8fH61atUp///vf5ePj45rfqjj35xu7w+Eo1KYqzf25goKC1KpVKyUkJFSL/+bPdvbYi1Ie804oq0DHjh1TUlKS6tevX9lduSDGGD322GP67LPPtHz5ckVHR7t9Hh0drXr16ik+Pt5VlpeXp1WrVqlTp04V3d1yd77xF6WqzH1RjDHKzc2t8vNelIKxF6Uqzfn111+vHTt2aPv27a5Xhw4ddN9992n79u1q0qRJlZ37843d29u7UJuqNPfnys3N1e7du1W/fv1q99/82WMvSrnMe5nPscFkZGSYbdu2mW3bthlJZubMmWbbtm1m//79JiMjwzz55JNm/fr1Zu/evWbFihUmNjbWNGzY0KSnp1d21y/II488YsLCwszKlStNcnKy65WVleWq8+KLL5qwsDDz2WefmR07dph77rnH1K9f/08/dmPOP/6qPPcTJkwwq1evNnv37jU///yzeeaZZ4yXl5f57rvvjDFVe95LGntVnvPinHsJryrP/bnOHntVn/snn3zS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",
      "text/plain": [
       "<Figure size 700x500 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.figure(figsize=(7,5))\n",
    "plt.scatter(patients[\"CAG_Repeats\"], patients[\"Disease\"], alpha=0.6)\n",
    "plt.xlabel(\"CAG Repeat\")\n",
    "plt.ylabel(\"Disease Label\")\n",
    "plt.title(\"Disease Prediction Pattern\")\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 78,
   "id": "63669236-9241-41d9-8fe3-9b21e44a3cb7",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAA7YAAAE6CAYAAAAiF84SAAAAOXRFWHRTb2Z0d2FyZQBNYXRwbG90bGliIHZlcnNpb24zLjkuMiwgaHR0cHM6Ly9tYXRwbG90bGliLm9yZy8hTgPZAAAACXBIWXMAAA9hAAAPYQGoP6dpAAAviUlEQVR4nO3deXhN5/7//9fOPEeEiBBinmlFqaENRakIHdWstFotRUc6nENbRbXn1GkPVa2hLS16OL5qjimtI61ZKUqVogQlISgy3L8//LI/tkSEkuSW5+O69nXZ93rvte617ntHXllrr+0wxhgBAAAAAGApt4LuAAAAAAAAfwXBFgAAAABgNYItAAAAAMBqBFsAAAAAgNUItgAAAAAAqxFsAQAAAABWI9gCAAAAAKxGsAUAAAAAWI1gCwAAAACwGsEWAPLgxx9/VO/evVWhQgX5+PgoICBA9evX15gxY3TixIkcX1O/fn05HA699957ua57/vz56tixoyIiIuTl5aXAwEDdfvvtGjZsmPbv33/Vvg0fPlwOh8P58PT0VLly5dS3b18lJSVd1/4WpO3bt2v48OHat2/fVWsv3e/cHqtWrbrp/c4ahz/++OOmb8tmjz32mAICAq64PCAgQI899pgkqXnz5nka37zWDR8+/Kr9+6vvRwBAwfAo6A4AQGH3ySef6JlnnlG1atX00ksvqWbNmkpLS9P69es1YcIEJSYm6r///a/LazZv3qxNmzZJkiZNmqQXX3wx23ozMzPVu3dvff7557rvvvs0atQoRUVF6c8//9S6des0ZcoUTZ48WQcOHMhTPxcvXqzg4GCdPn1aS5cu1T/+8Q+tWbNGmzdvlqen518/EPlk+/bteuONN9S8eXNFRUXlWpuYmOjy/K233tLKlSu1YsUKl/aaNWve6G4iH4wfP16nTp1yPl+wYIFGjBihKVOmqHr16s72CxcuyMvL66p1ZcuWveK2bvT7EQCQvwi2AJCLxMREPf3002rdurXmzp0rb29v57LWrVvrhRde0OLFi7O97tNPP5UkxcbGasGCBVqzZo2aNGniUvPOO+/o888/16hRozR06FCXZW3bttUrr7yijz/+OM99jY6OVokSJSRJrVq10h9//KEpU6Zo9erVatGiRZ7XY5M777zT5XnJkiXl5uaWrf1yZ8+elZ+f383sWqGSkZGh9PR0l/lrg8v/ILFz505JUu3atdWgQYMrvi6vdZe60e9HAED+4lJkAMjFyJEj5XA4NHHixBxDgZeXlzp06ODSdu7cOX355ZeKjo7W+++/L0maPHmyS82FCxc0ZswY1a5dO9sv0Vk8PDzUv3//6+571i/0R44ccWlftmyZWrZsqaCgIPn5+alp06Zavny5S03WZbWbNm3Sgw8+qKCgIAUHB6t79+46duxYtm3NnDlTjRs3lr+/vwICAtSmTRvnGess69evV+fOnRUVFSVfX19FRUWpS5cu+u2335w1U6dO1SOPPCJJatGihfMS0qlTp173cWjevLlq166tb7/9Vk2aNJGfn5/69Onj7Pe9996r0qVLy9fXVzVq1NDQoUN15syZbOv54YcfFBcXp9DQUPn4+KhSpUoaPHhwrtveuXOnKlasqEaNGuno0aNXrPvll1/Uu3dvValSRX5+fipTpozi4uK0devWbLUpKSl64YUXVLFiRXl7eyssLEzt2rVzhrl9+/bJ4XBozJgxGjFihCpUqCBvb2+tXLlSkjRv3jw1btxYfn5+CgwMVOvWrbOd+T527JiefPJJRUZGytvbWyVLllTTpk21bNkyZ82mTZvUvn17hYWFydvbWxEREYqNjdXBgwdzPSaF0V95P+Zl7mddfv3LL7+oXbt2CggIUGRkpF544QWdP38+W19GjBih6tWrO4997969c3zfAQD+D8EWAK4gIyNDK1asUHR0tCIjI/P8ujlz5ig5OVl9+vRRlSpV1KxZM82cOVOnT5921qxfv14pKSmKi4u7GV2XJO3du1eSVLVqVWfbtGnTdO+99yooKEifffaZZs2apeLFi6tNmzbZwq0kPfDAA6pcubL+85//aPjw4Zo7d67atGmjtLQ0Z83IkSPVpUsX1axZU7NmzdIXX3yh1NRU3XXXXdq+fbuzbt++fapWrZrGjh2rJUuW6J133tHhw4d1xx13OD+XGhsbq5EjR0qSxo0bp8TERCUmJio2NvYvHYvDhw+re/fu6tq1qxYuXKhnnnlGkrR79261a9dOkyZN0uLFizV48GDNmjUr27gsWbJEd911l/bv369//vOfWrRokV5//fVsfzS4VEJCgpo0aaK6detq5cqVCgsLu2LtoUOHFBoaqtGjR2vx4sUaN26cPDw81KhRI/3888/OutTUVDVr1kwff/yxevfurW+++UYTJkxQ1apVdfjwYZd1fvDBB1qxYoXee+89LVq0SNWrV9eXX36pjh07KigoSF999ZUmTZqk5ORkNW/eXKtXr3a+tkePHpo7d67+/ve/a+nSpfr000/VqlUrHT9+XJJ05swZtW7dWkeOHNG4ceMUHx+vsWPHqly5ckpNTc3TmKSnp+f4KAjX+37M69yXpLS0NHXo0EEtW7bU//t//099+vTR+++/r3feecdZk5mZqY4dO2r06NHq2rWrFixYoNGjRys+Pl7NmzfXn3/+eUP2FwBuSQYAkKOkpCQjyXTu3PmaXnfPPfcYHx8fk5ycbIwxZsqUKUaSmTRpkrNmxowZRpKZMGFCttenpaW5PK5m2LBhRpJJSkoyaWlpJjk52cyaNcv4+/ubLl26OOvOnDljihcvbuLi4lxen5GRYerVq2caNmyYbZ3PPfecS+306dONJDNt2jRjjDH79+83Hh4e5tlnn3WpS01NNeHh4aZTp05X7Hd6ero5ffq08ff3N//617+c7V9//bWRZFauXHnVfb9cr169jL+/v0tbTEyMkWSWL1+e62szMzNNWlqaSUhIMJLMli1bnMsqVapkKlWqZP78888rvj7rmB07dsx88cUXxsvLywwcONBkZGRc836kp6ebCxcumCpVqriMwZtvvmkkmfj4+Cu+du/evUaSqVSpkrlw4YKzPSMjw0RERJg6deq49Ck1NdWEhYWZJk2aONsCAgLM4MGDr7iN9evXG0lm7ty517xvvXr1MpJyffTq1SvH12a9l9atW5frNvJal+V63o/XMvez9nnWrFkute3atTPVqlVzPv/qq6+MJDN79myXunXr1hlJZvz48XnaHwAoijhjCwA30N69e7Vy5Uo9+OCDKlasmCTpkUceUWBgYLbLkXOSkpIiT09Pl8f69evztO3w8HB5enoqJCREnTp1UnR0tD777DPn8jVr1ujEiRPq1auXyxmyzMxMtW3bVuvWrct2CW63bt1cnnfq1EkeHh7Oy1qXLFmi9PR09ezZ02WdPj4+iomJcbkb8enTpzVkyBBVrlxZHh4e8vDwUEBAgM6cOaMdO3bkaR+vV0hIiO65555s7b/++qu6du2q8PBwubu7y9PTUzExMZLk7NOuXbu0Z88ePf744/Lx8bnqtt5++2099thjGj16tP71r3/Jze3q/9Wmp6dr5MiRqlmzpry8vOTh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4sQJU758efPYY4+ZH374wezdu9csW7bM/PLLL84a5sCtbcSIESY0NNTMnz/f7N2713z99dcmICDAjB071lnDHLi1LFy40Lz22mtm9uzZRpL573//67I8v8Y7PT3d1K5d27Ro0cJs3LjRxMfHm4iICDNgwICbfgyKstzGPyUlxbRq1crMnDnT7Ny50yQmJppGjRqZ6Ohol3Uw/rgehTLYNmzY0PTr18+lrXr16mbo0KEF1CPcKEePHjWSTEJCgjHGmMzMTBMeHm5Gjx7trDl37pwJDg42EyZMMMZc/CHo6elpZsyY4az5/fffjZubm1m8eLExxpjt27cbSeb777931iQmJhpJZufOnfmxa8hFamqqqVKliomPjzcxMTHOYMv43/qGDBlimjVrdsXlzIFbX2xsrOnTp49L24MPPmi6d+9ujGEO3OouDzb5Od4LFy40bm5u5vfff3fWfPXVV8bb29ucPHnypuwvXOX0h43LrV271khynsBi/HG9Ct2lyBcuXNCGDRt07733urTfe++9WrNmTQH1CjfKyZMnJUnFixeXJO3du1dJSUku4+3t7a2YmBjneG/YsEFpaWkuNREREapdu7azJjExUcHBwWrUqJGz5s4771RwcDDzphDo37+/YmNj1apVK5d2xv/WN2/ePDVo0ECPPPKIwsLCdPvtt+uTTz5xLmcO3PqaNWum5cuXa9euXZKkLVu2aPXq1WrXrp0k5kBRk5/jnZiYqNq1aysiIsJZ06ZNG50/f97l4xAoWCdPnpTD4VCxYsUkMf64fh4F3YHL/fHHH8rIyFCpUqVc2kuVKqWkpKQC6hVuBGOMnn/+eTVr1ky1a9eWJOeY5jTev/32m7PGy8tLISEh2WqyXp+UlKSwsLBs2wwLC2PeFLAZM2Zo48aNWrduXbZljP+t79dff9VHH32k559/Xq+++qrWrl2rgQMHytvbWz179mQOFAFDhgzRyZMnVb16dbm7uysjI0Nvv/22unTpIomfA0VNfo53UlJStu2EhITIy8uLOVFInDt3TkOHDlXXrl0VFBQkifHH9St0wTaLw+FweW6MydYGuwwYMEA//vijVq9enW3Z9Yz35TU51TNvCtaBAwc0aNAgLV26VD4+PlesY/xvXZmZmWrQoIFGjhwpSbr99tv1008/6aOPPlLPnj2ddcyBW9fMmTM1bdo0ffnll6pVq5Y2b96swYMHKyIiQr169XLWMQeKlvwab+ZE4ZWWlqbOnTsrMzNT48ePv2o944+rKXSXIpcoUULu7u7Z/pJy9OjRbH91gT2effZZzZs3TytXrlTZsmWd7eHh4ZKU63iHh4frwoULSk5OzrXmyJEj2bZ77Ngx5k0B2rBhg44eParo6Gh5eHjIw8NDCQkJ+uCDD+Th4eEcG8b/1lW6dGnVrFnTpa1GjRrav3+/JH4GFAUvvfSShg4dqs6dO6tOnTrq0aOHnnvuOY0aNUoSc6Coyc/xDg8Pz7ad5ORkpaWlMScKWFpamjp16qS9e/cqPj7eebZWYvxx/QpdsPXy8lJ0dLTi4+Nd2uPj49WkSZMC6hWulzFGAwYM0Jw5c7RixQpVqFDBZXmFChUUHh7uMt4XLlxQQkKCc7yjo6Pl6enpUnP48GFt27bNWdO4cWOdPHlSa9euddb88MMPOnnyJPOmALVs2VJbt27V5s2bnY8GDRqoW7du2rx5sypWrMj43+KaNm2a7Su+du3apfLly0viZ0BRcPbsWbm5uf664e7u7vy6H+ZA0ZKf4924cWNt27ZNhw8fdtYsXbpU3t7eio6Ovqn7iSvLCrW7d+/WsmXLFBoa6rKc8cd1y887VeVV1tf9TJo0yWzfvt0MHjzY+Pv7m3379hV013CNnn76aRMcHGxWrVplDh8+7HycPXvWWTN69GgTHBxs5syZY7Zu3Wq6dOmS423/y5Yta5YtW2Y2btxo7rnnnhxv+163bl2TmJhoEhMTTZ06dfiah0Lo0rsiG8P43+rWrl1rPDw8zNtvv212795tpk+fbvz8/My0adOcNcyBW1uvXr1MmTJlnF/3M2fOHFOiRAnz8ssvO2uYA7eW1NRUs2nTJrNp0yYjyfzzn/80mzZtct71Nr/GO+vrXlq2bGk2btxoli1bZsqWLcvXvdxkuY1/Wlqa6dChgylbtqzZvHmzy++G58+fd66D8cf1KJTB1hhjxo0bZ8qXL2+8vLxM/fr1nV8PA7tIyvExZcoUZ01mZqYZNmyYCQ8PN97e3ubuu+82W7dudVnPn3/+aQYMGGCKFy9ufH19Tfv27c3+/ftdao4fP266detmAgMDTWBgoOnWrZtJTk7Oh73Etbg82DL+t75vvvnG1K5d23h7e5vq1aubiRMnuixnDtzaTp06ZQYNGmTKlStnfHx8TMWKFc1rr73m8kssc+DWsnLlyhz/7+/Vq5cxJn/H+7fffjOxsbHG19fXFC9e3AwYMMCcO3fuZu5+kZfb+O/du/eKvxuuXLnSuQ7GH9fDYYwx+Xd+GAAAAACAG6vQfcYWAAAAAIBrQbAFAAAAAFiNYAsAAAAAsBrBFgAAAABgNYItAAAAAMBqBFsAAAAAgNUItgAAAAAAqxFsAQAAAABWI9gCAPAXNG/eXIMHDy7obuRo1apVcjgcSklJybUuKipKY8eOzZc+AQBwMxBsAQAFJikpSYMGDVLlypXl4+OjUqVKqVmzZpowYYLOnj1b0N3Lkzlz5uitt9667tc3b95cDodDDodD3t7eqlq1qkaOHKmMjIy/3LcmTZro8OHDCg4OliRNnTpVxYoVy1a3bt06Pfnkk395ewAAFBSPgu4AAKBo+vXXX9W0aVMVK1ZMI0eOVJ06dZSenq5du3Zp8uTJioiIUIcOHQq6m1dVvHjxv7yOvn376s0339S5c+c0f/58DRw4UO7u7hoyZMhfWq+Xl5fCw8OvWleyZMm/tB0AAAoaZ2wBAAXimWeekYeHh9avX69OnTqpRo0aqlOnjh566CEtWLBAcXFxztqTJ0/qySefVFhYmIKCgnTPPfdoy5YtzuXDhw/Xbbfdpi+++EJRUVEKDg5W586dlZqa6qw5f/68Bg4cqLCwMPn4+KhZs2Zat26dc3nWZbtLlizR7bffLl9fX91zzz06evSoFi1apBo1aigoKEhdunRxOZt8+aXI58+f18svv6zIyEh5e3urSpUqmjRpUq7Hws/PT+Hh4YqKitKAAQPUsmVLzZ07V5KUnJysnj17KiQkRH5+frrvvvu0e/du52t/++03xcXFKSQkRP7+/qpVq5YWLlzosk8pKSlatWqVevfurZMnTzrPEA8fPlxS9kuR9+/fr44dOyogIEBBQUHq1KmTjhw5ck3HGwCA/ESwBQDku+PHj2vp0qXq37+//P39c6xxOBySJGOMYmNjlZSUpIULF2rDhg2qX7++WrZsqRMnTjjr9+zZo7lz52r+/PmaP3++EhISNHr0aOfyl19+WbNnz9Znn32mjRs3qnLlymrTpo3LOqSLoe3f//631qxZowMHDqhTp04aO3asvvzySy1YsEDx8fH68MMPr7hvPXv21IwZM/TBBx9ox44dmjBhggICAq7p+Pj6+iotLU2S9Nhjj2n9+vWaN2+eEhMTZYxRu3btnMv79++v8+fP69tvv9XWrVv1zjvv5Li9Jk2aaOzYsQoKCtLhw4d1+PBhvfjii9nqjDG6//77deLECSUkJCg+Pl579uzRo48+6lJ3teMNAEC+MgAA5LPvv//eSDJz5sxxaQ8NDTX+/v7G39/fvPzyy8YYY5YvX26CgoLMuXPnXGorVapkPv74Y2OMMcOGDTN+fn7m1KlTzuUvvfSSadSokTHGmNOnTxtPT08zffp05/ILFy6YiIgIM2bMGGOMMStXrjSSzLJly5w1o0aNMpLMnj17nG1PPfWUadOmjfN5TEyMGTRokDHGmJ9//tlIMvHx8Xk+Fpe+PiMjwyxatMh4eXmZl19+2ezatctIMv/73/+c9X/88Yfx9fU1s2bNMsYYU6dOHTN8+PAc1521T8nJycYYY6ZMmWKCg4Oz1ZUvX968//77xhhjli5datzd3c3+/fudy3/66Scjyaxdu9YYc/XjDQBAfuOMLQCgwGSdlc2ydu1abd68WbVq1dL58+clSRs2bNDp06cVGhqqgIAA52Pv3r3as2eP87VRUVEKDAx0Pi9durSOHj0q6eLZxbS0NDVt2tS53NPTUw0bNtSOHTtc+lC3bl3nv0uVKiU/Pz9VrFjRpS1rvZfbvHmz3N3dFRMTc03HYfz48QoICJCPj486dOig7t27a9iwYdqxY4c8PDzUqFEjZ21oaKiqVavm7PfAgQM1YsQINW3aVMOGDdOPP/54Tdu+3I4dOxQZGanIyEhnW82aNVWsWDGXY5Xb8QYAIL9x8ygAQL6rXLmyHA6Hdu7c6dKeFSB9fX2dbZmZmSpdurRWrVqVbT2X3uHX09PTZZnD4VBmZqaki5fXZrVdyhiTre3S9TgcjlzXe7lL+30tunXrptdee03e3t6KiIiQu7u7S78vd2m/n3jiCbVp00YLFizQ0qVLNWrUKP3jH//Qs88+e119yemY5NR+LccFAICbjTO2AIB8FxoaqtatW+vf//63zpw5k2tt/fr1lZSUJA8PD1WuXNnlUaJEiTxtr3LlyvLy8tLq1audbWlpaVq/fr1q1Kjxl/blUnXq1FFmZqYSEhKu6XXBwcGqXLmyIiMjnaFWunimND09XT/88IOz7fjx49q1a5dLvyMjI9WvXz/NmTNHL7zwgj755JMct+Pl5XXVrxGqWbOm9u/frwMHDjjbtm/frpMnT97QYwUAwI1EsAUAFIjx48crPT1dDRo00MyZM7Vjxw79/PPPmjZtmnbu3OkMeK1atVLjxo11//33a8mSJdq3b5/WrFmj119/XevXr8/Ttvz9/fX000/rpZde0uLFi7V9+3b17dtXZ8+e1eOPP37D9ikqKkq9evVSnz59NHfuXO3du1erVq3SrFmzrmt9VapUUceOHdW3b1+tXr1aW7ZsUffu3VWmTBl17NhRkjR48GAtWbJEe/fu1caNG7VixYorBtCoqCidPn1ay5cv1x9//JHjdwW3atVKdevWVbdu3bRx40atXbtWPXv2VExMjBo0aHBd+wEAwM1GsAUAFIhKlSpp06ZNatWqlV555RXVq1dPDRo00IcffqgXX3xRb731lqSLl7guXLhQd999t/r06aOqVauqc+fO2rdvn0qVKpXn7Y0ePVoPPfSQevToofr16+uXX37RkiVLFBISckP366OPPtLDDz+sZ555RtWrV1ffvn2velY6N1OmTFF0dLTat2+vxo0byxijhQsXOi8FzsjIUP/+/VWjRg21bdtW1apV0/jx43NcV5MmTdSvXz89+uijKlmypMaMGZOtxuFwaO7cuQoJCdHdd9+tVq1aqWLFipo5c+Z17wMAADebw1zpAzwAAAAAAFiAM7YAAAAAAKsRbAEAAAAAViPYAgAAAACsRrAFAAAAAFiNYAsAAAAAsBrBFgAAAABgNYItAAAAAMBqBFsAAAAAgNUItgAAAAAAqxFsAQAAAABWI9gCAAAAAKz2/wEQAdR+TvbZNwAAAABJRU5ErkJggg==",
      "text/plain": [
       "<Figure size 1200x300 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.figure(figsize=(12,3))\n",
    "\n",
    "for start, end, length in positions:\n",
    "    plt.plot([start, end], [1,1], linewidth=length, color=\"blue\")\n",
    "\n",
    "plt.title(\"CAG Repeat Track across HTT Gene\")\n",
    "plt.xlabel(\"Genomic Position\")\n",
    "plt.yticks([])\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 80,
   "id": "eed1a8c3-68bb-47ea-84ba-0a3e61ac2826",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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",
      "text/plain": [
       "<Figure size 300x800 with 2 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "import seaborn as sns\n",
    "\n",
    "heat = np.array(repeat_counts).reshape(-1,1)\n",
    "\n",
    "plt.figure(figsize=(3,8))\n",
    "sns.heatmap(heat, cmap=\"Reds\", cbar=True)\n",
    "plt.title(\"CAG Repeat Length Heatmap\")\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 82,
   "id": "cdfe824e-7472-49c7-b55c-e8ab82490310",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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PzBcad26dTp27JgSExPtRw0lFenz5nOTfeO31atX5zrf29u7UNcXExOj3bt366677tI//vEPtW/fXtWrV8/RLrc6HT9+3B6kzp49qxUrVigmJkYjR460t7l48aL+93//94bH9+CDD+rBBx/U5cuXtWXLFs2cOVNDhgxRYGCgunXrlusyfn5+Sk5OzjE9+0Z6hXk3f09PT3Xv3l3/+te/lJycrLlz58rb21t/+9vfrrtsgwYNtGTJElmWpV27dmnevHkaN26cPD09NXLkyAL93S9cuFB///vfNX78eIfpf/755w0/hvG5557T119/rVWrVqlRo0YO827V5+l21g4ASgKOyANAKda6dWtJf4WBK3322WdKS0uzz7+WF154Qbt27bLf4OvZZ591mB8eHq49e/bkOIV4yZIl1+07MDBQHh4e2rVrl8P0zz//PEfbvI7IZZ8Wnj0/27vvvnvd9d9O7du316lTp3T58uUcR9SbNWumOnXqFNq61qxZowkTJuiVV17RmjVr5OvrqyeffDLXG8ktXrxYlmXZ3x8+fFjr169XRESEpL/2r2VZOfbvnDlzrnvTuPxwdnbWPffcY78D+rZt2/Js27p1a3vQvNKCBQvk5eVV6I+re/rpp3X58mVNnjxZq1atUrdu3eTl5ZXv5W02mxo1aqQ333xT5cqVs29bQf7ubTZbjn2/cuVK/f777zewRdIrr7yiuLg4zZkzRw8//HCu65Py93kqyFHy1q1ba+/evTnqu2DBAtlsNrVs2TLf2wAApQFH5AGgFHvkkUcUGRmpESNGKCUlRffff7927dqlmJgYNWnSRL169cpXH2FhYUpISLA/OupKQ4YM0dy5c9W2bVuNGzdOgYGBWrRokfbv3y9JcnLK+zdlm82mnj17au7cuapRo4YaNWqkTZs2adGiRTnaZt9tffr06erdu7dcXV1Vp04dtWjRQuXLl1e/fv0UExMjV1dXffjhh9q5c2dBdtUt161bN3344YeKiorSCy+8oLvvvluurq46evSoEhIS1LFjR3Xu3Pm6/WzdulW+vr45poeFhcnHx0fJycnq2bOnwsPDFRMTIycnJ3300Ud66KGH9PLLL2vatGkOy508eVKdO3fWs88+q7NnzyomJkYeHh4aNWqUJMnHx0cPPfSQJk+eLH9/f4WGhiopKUnvv//+DR8Rfuedd7Ru3Tq1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GIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAAABiEIA8AAAAAgEEI8gAAAAAAGIQgDwAAAACAQQjyAAAAAAAYhCAPAAAAAIBBCPIAAAAA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",
      "text/plain": [
       "<Figure size 1200x400 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.figure(figsize=(12,4))\n",
    "\n",
    "plt.plot([196,259],[2,2], linewidth=normal_max, color=\"green\", label=\"Normal\")\n",
    "\n",
    "plt.plot([196,196+45*3],[1,1], linewidth=mutant_max, color=\"red\", label=\"Mutant\")\n",
    "\n",
    "plt.title(\"Polyglutamine Expansion Visualization\")\n",
    "plt.yticks([1,2],[\"Mutant\",\"Normal\"])\n",
    "plt.legend()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 84,
   "id": "498ab2b5-433e-466c-8a8a-acbe01058941",
   "metadata": {},
   "outputs": [],
   "source": [
    "from sklearn.metrics import roc_curve, auc\n",
    "\n",
    "y_prob = model.predict_proba(X_test)[:,1]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 86,
   "id": "38541ed7-7286-4b3b-8c0b-0d940cd78450",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "1.0"
      ]
     },
     "execution_count": 86,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "fpr, tpr, thresholds = roc_curve(y_test, y_prob)\n",
    "\n",
    "roc_auc = auc(fpr, tpr)\n",
    "\n",
    "roc_auc"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 88,
   "id": "a0d62ab6-f441-4c1b-a3e4-278e12db0f85",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAhgAAAHUCAYAAAB4RlFCAAAAOXRFWHRTb2Z0d2FyZQBNYXRwbG90bGliIHZlcnNpb24zLjkuMiwgaHR0cHM6Ly9tYXRwbG90bGliLm9yZy8hTgPZAAAACXBIWXMAAA9hAAAPYQGoP6dpAAByAUlEQVR4nO3dd1QU198G8GdhYekIIiCKiAV7BQsSe0FQwI69ocZoYtdoTGJJIUZjjIli12gUMHaNJVhibyB2jUZRLKCCCkhfuO8fvuzPFVAWF4byfM7Zc9i7U54Zdtkvd+7MyIQQAkRERERapCN1ACIiIip5WGAQERGR1rHAICIiIq1jgUFERERaxwKDiIiItI4FBhEREWkdCwwiIiLSOhYYREREpHUsMIiIiEjrWGAUU+vWrYNMJlM95HI5ypcvj759++L27ds5zpOeno6AgAC4urrC3NwchoaGqFWrFqZPn47Y2Ngc58nMzMSGDRvQoUMHWFlZQU9PD9bW1ujatSt2796NzMzM92ZNTU3Fb7/9ho8++ggWFhbQ19dHhQoV0KdPHxw9evSD9kNJM3v2bMhkMsTExOT4et26ddGmTZsCzZCUlITZs2fjn3/+yfZa1vvu3r17BZph06ZNWLRoUYGu413u3bun9vnS09ND2bJl0aRJE0ycOBHXrl3LNs8///wDmUyW434rad7cNzKZDObm5mjTpg3++uuvQll/1ufkTZUrV8bQoUM1Wk5ReK+XZCwwirm1a9fi9OnTOHjwID799FPs2rULH330EV68eKE2XVJSEjp27IjPPvsMjRo1QmBgIPbu3YtBgwZhxYoVaNSoEf7991+1eVJSUuDp6YkhQ4bA2toaAQEBOHz4MJYtWwY7Ozv07t0bu3fvfme+mJgYuLm5YdKkSahbty7WrVuHQ4cO4aeffoKuri7at2+PS5cuaX2/UP4lJSVhzpw5Of7R7dKlC06fPo3y5csXaAapC4wsn332GU6fPo2jR49iw4YN6NatG3bt2oUGDRpg/vz5atM2btwYp0+fRuPGjSVKW7h69eqF06dP4+TJk1iyZAmio6Ph5eVVaEXG27Zv346vvvpKo3mKwnu9RBNULK1du1YAEOfPn1drnzNnjgAg1qxZo9Y+atQoAUAEBQVlW9a///4rzM3NRZ06dYRSqVS1f/LJJwKA+P3333PMcOvWLXHp0qV35vTw8BByuVwcOnQox9fPnTsn7t+//85l5FVSUpJWliOlWbNmCQDi2bNnOb5ep04d0bp16wLN8OzZMwFAzJo1q0DX8y5dunQRDg4Okq0/IiJCABDz58/P9lpSUpLo3LmzACD27t0rQTrpARBjx45Va/vvv/8EANGhQ4dc50tLSxPp6ekfvP6sz8mHKgrv9ZKMPRgljIuLCwDgyZMnqrbo6GisWbMG7u7u8PX1zTaPk5MTPv/8c1y7dg07duxQzbNq1Sq4u7tj8ODBOa6revXqqF+/fq5ZwsLCsG/fPvj5+aFdu3Y5TtOkSRNUqlQJQM7dnkDOXZWVK1dG165dsW3bNjRq1AgGBgaYM2cOGjVqhJYtW2ZbRkZGBipUqIAePXqo2tLS0vDtt9+iZs2aUCgUKFeuHIYNG4Znz57luk1FTW7d8lld/OvWrVO1DR06FCYmJvjvv//g6ekJExMT2NvbY/LkyUhNTVXNV65cOQDAnDlzVF3gWV3POf0u2rRpg7p16+L8+fNo2bIljIyMUKVKFfzwww/ZDqFdu3YNnTp1gpGREcqVK4exY8fir7/+UtuGrK72+/fvq3XDZ3n+/DnGjBmDChUqQF9fH1WqVMHMmTNV25BFJpPh008/xYYNG1CrVi0YGRmhQYMG2LNnzwfsccDQ0BCrV6+Gnp6eWi9GTr+Lu3fvom/fvrCzs4NCoYCNjQ3at2+Pixcvqi0zODgYrq6uMDY2homJCdzd3REeHq42TWhoKPr27YvKlSvD0NAQlStXRr9+/XD//n216ZKSkjBlyhQ4OjrCwMAAlpaWcHFxQWBgYLbleXt7w9LSEgYGBmjUqBE2b96c7/1StWpVlCtXTpUna39s2LABkydPRoUKFaBQKPDff/8BAA4ePIj27dvDzMwMRkZGcHNzw6FDh7It96+//kLDhg2hUCjg6OiIBQsW5Lj+nA6RvHz5EpMnT0aVKlWgUChgbW0NT09P3Lx5M1/vdQBYs2YNGjRooNq33bt3x40bN9SmyctnrTSQSx2AtCsiIgLA66Ihy5EjR6BUKtGtW7dc5+vWrRu++OILhISEoGfPnjhy5AjS09PfOc/7/P3336plF4QLFy7gxo0b+PLLL+Ho6AhjY2PY2dlh/PjxuH37NqpXr66W5fHjxxg2bBiA12NLfHx8cPz4cUybNg0tWrTA/fv3MWvWLLRp0wahoaEwNDQskNx5kZGRAaVSqfXlpqenw9vbG35+fpg8eTKOHTuGb775Bubm5vj6669Rvnx57N+/H507d4afnx9GjBgBAKo/xLmJjo7GgAEDMHnyZMyaNQvbt2/HjBkzYGdnpypQo6Ki0Lp1axgbGyMgIADW1tYIDAzEp59+qraspUuXYtSoUbhz5w62b9+u9lpKSgratm2LO3fuYM6cOahfvz6OHz8Of39/XLx4MVv3/F9//YXz589j7ty5MDExwY8//oju3bvj33//RZUqVfK9H+3s7ODs7IxTp05BqVRCLs/5T6mnpycyMjLw448/olKlSoiJicGpU6fw8uVL1TTff/89vvzySwwbNgxffvkl0tLSMH/+fLRs2RLnzp1D7dq1Abwu/mrUqIG+ffvC0tISUVFRCAgIQJMmTXD9+nVYWVkBACZNmoQNGzbg22+/RaNGjZCYmIirV6+qjbM6cuQIOnfujGbNmmHZsmUwNzdHUFAQfH19kZSUpPFYBgB48eIFYmNj1T53ADBjxgy4urpi2bJl0NHRgbW1Nf744w8MHjwYPj4++P3336Gnp4fly5fD3d0dBw4cQPv27QEAhw4dgo+PD1xdXREUFKTal2/+A5WbhIQEfPTRR7h37x4+//xzNGvWDK9evcKxY8cQFRWFFi1aaPxe9/f3xxdffIF+/frB398fsbGxmD17NlxdXXH+/Hm1bX/fZ61UkLoLhfIn6xDJmTNnRHp6ukhISBD79+8Xtra2olWrVmrdkD/88IMAIPbv35/r8pKTkwUA4eHhked53mf06NECgLh582aeps+t2zNrWyMiIlRtDg4OQldXV/z7779q08bExAh9fX3xxRdfqLX36dNH2NjYqPZLYGCgACC2bt2qNt358+cFALF06dI8Zda2rH3wrsebh0iOHDkiAIgjR46oLSeri3/t2rWqtiFDhggAYvPmzWrTenp6iho1aqiev6vbOKffRevWrQUAcfbsWbVpa9euLdzd3VXPp06dKmQymbh27ZradO7u7tm2IbdDJMuWLctxG+bNmycAiL///lvVBkDY2NiI+Ph4VVt0dLTQ0dER/v7+2Zb9pncdIsni6+srAIgnT54IIbL/LmJiYgQAsWjRolyXERkZKeRyufjss8/U2hMSEoStra3o06dPrvMqlUrx6tUrYWxsLH755RdVe926dUW3bt3euX01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3L59G7t27cJnn32W74xTp07FunXrsGzZMty+fRsLFy7Etm3bVIMb161bh9WrV+Pq1auqbTA0NISDg0OOy6tcuTIOHTqE6OhovHjxItf1DhgwABcuXMB3332HXr16wcDAQPWatrbTzMwMI0aMwKxZsyCEyNPvo3Llynj16hUOHTqEmJgYJCUlwcnJCQMGDMDgwYOxbds2RERE4Pz585g3bx727t2rUSaiIkvKASBElLMhQ4YIHx+fHF9buHChKF++vDA0NBTu7u5i/fr1AoB48eKFEEJ9UGFqaqro27evsLe3F/r6+sLOzk58+umnagMbz507Jzp27ChMTEyEsbGxqF+/vvjuu+9yzZbToMW3LV26VFSpUkXo6ekJJycnsX79etVr27dvF82aNRNmZmbC2NhYNG/eXBw8eFD1+tuDPHft2iWqVasm5HK5cHBwEEJkH+SZpUmTJgKAOHz4cLbXtLWd9+/fF3K5XAQHBwsh3v/7EEKI0aNHi7JlywoAYtasWUIIIdLS0sTXX38tKleuLPT09IStra3o3r27uHz5cq6ZiIoTmRBCSFviEBERUUnDQyRERESkdSwwiIiISOtYYBAREZHWscAgIiIirWOBQURERFrHAoOIiIi0jgUGERERaR0LDCIiItI6FhhERESkdSwwiIiISOtYYBAREZHW/R9SCbAzDVr8rQAAAABJRU5ErkJggg==",
      "text/plain": [
       "<Figure size 600x500 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.figure(figsize=(6,5))\n",
    "\n",
    "plt.plot(fpr, tpr, color=\"darkorange\",\n",
    "         label=\"ROC curve (AUC = %0.2f)\" % roc_auc)\n",
    "\n",
    "plt.plot([0,1],[0,1], color=\"navy\", linestyle=\"--\")\n",
    "\n",
    "plt.xlabel(\"False Positive Rate\")\n",
    "plt.ylabel(\"True Positive Rate\")\n",
    "plt.title(\"ROC Curve — Huntington Disease Prediction\")\n",
    "plt.legend(loc=\"lower right\")\n",
    "\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 90,
   "id": "72c27ab7-220e-4905-ba73-339ccfa3dd27",
   "metadata": {},
   "outputs": [],
   "source": [
    "y_pred = model.predict(X_test)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 92,
   "id": "860abd7a-6571-447d-88cc-319153e3892f",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "array([[18,  0],\n",
       "       [ 0, 22]], dtype=int64)"
      ]
     },
     "execution_count": 92,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "from sklearn.metrics import confusion_matrix\n",
    "\n",
    "cm = confusion_matrix(y_test, y_pred)\n",
    "\n",
    "cm"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 94,
   "id": "0b42ecf6-91a9-442e-bb95-e3efa3b9daf6",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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keKKf+eKfp9Wr15d4qpZJycn1YpRbZ9BeT6XunXraizIq4jn/+YVZ+i/eXLxQk2BhoeHQwiBSZMmaV00kpeXh8OHDwMAevXqBQCqRSxFEhISkJycjN69exssrqKVjD/99JNaeVEs2piamsLX11e1Yu7SpUsl1u3duzdOnz6tSnhFdu7ciRo1alTa0vZ69eph7ty5CAoKQnBwcIn1ipaNP/+v65ycHHz++ecadQ01qi4oKMCoUaMgSRKOHTuGyMhIrF69Gvv37y+1nYWFBTp06KDXVp7pVV04OTnB0tJS43tz6NChCvVr6JkL4NmU+J9//oljx46ple/du1fr8QHNUWnR71zx38mYmBhkZ2cb9HeyIrp06YJatWrh119/LfG7YGFhAWtra3Tq1An79+9XGzX9+eefpf7uFwkMDMR///tf1dSwNrqM8P38/GBlZaXx+d69e1d1GoV088KMAIFnP+D169dj2rRpaN++Pd588020bNkSeXl5uHz5MjZu3IhWrVohKCgITZs2xeTJk7F69WqYmJggMDAQt27dwsKFC9GgQQO1lXMVNWDAANjb2yM0NBQffvghzMzMsH37do3prg0bNuD06dMYOHAg3NzckJubq1ppWdp5gIiICHz99dfo2bMnPvjgA9jb22P37t04cuQIli9fDjs7O4O9l+KWLl1aZp2BAwfik08+wejRozF58mQ8fPgQK1as0Dpiat26Nfbu3Yvo6Gh4enrC0tJSr8QSERGBb7/9FidOnICzszNmz56NuLg4hIaGwsfHR+9pzqokSRLGjh2LrVu3wsvLC97e3vjhhx/wxRdfVKjf1q1b4+zZszh8+DBcXFxgY2ODpk2bVqjP4OBgrFy5EmPHjsU///lPNGrUCMeOHcN//vMfAICJyd//Vi76eS5btgyBgYEwNTVFmzZt0KdPH/Tr1w/z589HVlYWunTpgp9++gkRERHw8fF5YS5fqVmzJlavXo3g4GBkZGRg2LBhcHR0xP3793HlyhXcv38f69evBwB89NFH6N+/v+o65IKCAixbtgzW1tZl3pVq5syZiI6OxqBBg/Duu++iU6dOyMnJQVxcHF599VX07NkTNjY2cHd3x6FDh9C7d2/Y29ujTp06Wu8oVKtWLSxcuBDvvfcexo8fj1GjRuHhw4dYvHgxLC0tERERURkfV/Vm7FU42iQlJYng4GDh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AM+mEn19fdX2FU3Vikp4vOaFCxcwcuRILF68GP369YOdnR327t2Lf//736W2kyRJ9b/liZ1IbpgAqcpZW1ujUaNG5a7frl07pKenw8zMDA0bNtRap3nz5rhw4QLGjx+vKrtw4UKJfTZu3BhWVlb45ptv8MYbb2jsLzrnV1BQoCpzcnJCvXr18Ntvv2HMmDFa+23RogU+//xz5OTkqJJsaXGUx/fffw93d3csWLBAVXb79m2Neqmpqfjjjz/g6uoKAIiPj4eJiQmaNGlSrtiJ5IYJkF54r7zyCvz8/DB48GAsW7YMTZs2xR9//IGjR49i8ODB6NChA2bMmIHg4GB06NABXbt2xe7du3H16tUSF8FYWlpi/vz5mDdvHiwsLNClSxfcv38fV69eRWhoKBwdHWFlZYXjx4+jfv36sLS0hJ2dHRYtWoR33nkHtra2CAwMhFKpRGJiIh49eoSwsDCMHj0aCxYsQGhoKN5//33cunULK1asKNf7vH//vsZ1h87OzmjUqBFSU1Oxd+9edOzYEUeOHMGBAwe0vqfg4GCsWLECWVlZeOeddzB8+HA4OzsDQJmxE8mOsU9CkrwUXwRTXEREhNrClSJZWVni7bffFq6ursLc3Fw0aNBAjBkzRqSmpqrqfPzxx6JOnTqiZs2aIjg4WMybN6/ERTBCCFFQUCD++c9/Cnd3d2Fubi7c3NzEkiVLVPs3bdokGjRoIExMTESPHj1U5bt37xZt27YVFhYWonbt2qJ79+5i//79qv3x8fHC29tbWFhYiLZt24qYmJhyLYIBoLFFREQIIYSYO3eucHBwEDVr1hQjRowQK1euFHZ2dhqf27p164Srq6uwtLQUQ4YMERkZGWrHKS12LoIhuZGEqISTFkRERC84XghPRESyxARIRESyxARIRESyxARIRESyxARIRESyxARIRESyxARIRESyxARIRESyxARIRESyxARIRESyxARIRESy9P8p0KzRZcqWigAAAABJRU5ErkJggg==",
      "text/plain": [
       "<Figure size 500x400 with 2 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "import seaborn as sns\n",
    "import matplotlib.pyplot as plt\n",
    "\n",
    "plt.figure(figsize=(5,4))\n",
    "sns.heatmap(cm, annot=True, fmt=\"d\", cmap=\"Blues\")\n",
    "\n",
    "plt.xlabel(\"Predicted Label\")\n",
    "plt.ylabel(\"True Label\")\n",
    "plt.title(\"Confusion Matrix — Huntington Prediction\")\n",
    "\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 96,
   "id": "11a16a3b-7791-444b-bbca-f2187af31101",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "              precision    recall  f1-score   support\n",
      "\n",
      "           0       1.00      1.00      1.00        18\n",
      "           1       1.00      1.00      1.00        22\n",
      "\n",
      "    accuracy                           1.00        40\n",
      "   macro avg       1.00      1.00      1.00        40\n",
      "weighted avg       1.00      1.00      1.00        40\n",
      "\n"
     ]
    }
   ],
   "source": [
    "from sklearn.metrics import classification_report\n",
    "\n",
    "print(classification_report(y_test, y_pred))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 98,
   "id": "2b3437d0-825c-40d5-9af4-918e56eb84f5",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "['huntington_ml_model.pkl']"
      ]
     },
     "execution_count": 98,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "import joblib\n",
    "\n",
    "joblib.dump(model, \"huntington_ml_model.pkl\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 100,
   "id": "5dfcd43f-cb96-446f-88e8-8ed160482da5",
   "metadata": {},
   "outputs": [],
   "source": [
    "loaded_model = joblib.load(\"huntington_ml_model.pkl\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 102,
   "id": "1421b143-4c36-4796-a220-566d369a8887",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "array([1], dtype=int64)"
      ]
     },
     "execution_count": 102,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "loaded_model.predict(pd.DataFrame({\"CAG_Repeats\":[42]}))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 104,
   "id": "465d5ed4-0ea5-4fbe-8e16-259bb7f1aa0b",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Requirement already satisfied: streamlit in c:\\programdata\\anaconda3\\lib\\site-packages (1.37.1)Note: you may need to restart the kernel to use updated packages.\n",
      "\n",
      "Requirement already satisfied: altair<6,>=4.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (5.0.1)\n",
      "Requirement already satisfied: blinker<2,>=1.0.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (1.6.2)\n",
      "Requirement already satisfied: cachetools<6,>=4.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (5.3.3)\n",
      "Requirement already satisfied: click<9,>=7.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (8.1.7)\n",
      "Requirement already satisfied: numpy<3,>=1.20 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (1.26.4)\n",
      "Requirement already satisfied: packaging<25,>=20 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (24.1)\n",
      "Requirement already satisfied: pandas<3,>=1.3.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (2.2.2)\n",
      "Requirement already satisfied: pillow<11,>=7.1.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (10.4.0)\n",
      "Requirement already satisfied: protobuf<6,>=3.20 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (4.25.3)\n",
      "Requirement already satisfied: pyarrow>=7.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (16.1.0)\n",
      "Requirement already satisfied: requests<3,>=2.27 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (2.32.3)\n",
      "Requirement already satisfied: rich<14,>=10.14.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (13.7.1)\n",
      "Requirement already satisfied: tenacity<9,>=8.1.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (8.2.3)\n",
      "Requirement already satisfied: toml<2,>=0.10.1 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (0.10.2)\n",
      "Requirement already satisfied: typing-extensions<5,>=4.3.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (4.11.0)\n",
      "Requirement already satisfied: gitpython!=3.1.19,<4,>=3.0.7 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (3.1.43)\n",
      "Requirement already satisfied: pydeck<1,>=0.8.0b4 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (0.8.0)\n",
      "Requirement already satisfied: tornado<7,>=6.0.3 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (6.4.1)\n",
      "Requirement already satisfied: watchdog<5,>=2.1.5 in c:\\programdata\\anaconda3\\lib\\site-packages (from streamlit) (4.0.1)\n",
      "Requirement already satisfied: jinja2 in c:\\programdata\\anaconda3\\lib\\site-packages (from altair<6,>=4.0->streamlit) (3.1.4)\n",
      "Requirement already satisfied: jsonschema>=3.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from altair<6,>=4.0->streamlit) (4.23.0)\n",
      "Requirement already satisfied: toolz in c:\\programdata\\anaconda3\\lib\\site-packages (from altair<6,>=4.0->streamlit) (0.12.0)\n",
      "Requirement already satisfied: colorama in c:\\programdata\\anaconda3\\lib\\site-packages (from click<9,>=7.0->streamlit) (0.4.6)\n",
      "Requirement already satisfied: gitdb<5,>=4.0.1 in c:\\programdata\\anaconda3\\lib\\site-packages (from gitpython!=3.1.19,<4,>=3.0.7->streamlit) (4.0.7)\n",
      "Requirement already satisfied: python-dateutil>=2.8.2 in c:\\programdata\\anaconda3\\lib\\site-packages (from pandas<3,>=1.3.0->streamlit) (2.9.0.post0)\n",
      "Requirement already satisfied: pytz>=2020.1 in c:\\programdata\\anaconda3\\lib\\site-packages (from pandas<3,>=1.3.0->streamlit) (2024.1)\n",
      "Requirement already satisfied: tzdata>=2022.7 in c:\\programdata\\anaconda3\\lib\\site-packages (from pandas<3,>=1.3.0->streamlit) (2023.3)\n",
      "Requirement already satisfied: charset-normalizer<4,>=2 in c:\\programdata\\anaconda3\\lib\\site-packages (from requests<3,>=2.27->streamlit) (3.3.2)\n",
      "Requirement already satisfied: idna<4,>=2.5 in c:\\programdata\\anaconda3\\lib\\site-packages (from requests<3,>=2.27->streamlit) (3.7)\n",
      "Requirement already satisfied: urllib3<3,>=1.21.1 in c:\\programdata\\anaconda3\\lib\\site-packages (from requests<3,>=2.27->streamlit) (2.2.3)\n",
      "Requirement already satisfied: certifi>=2017.4.17 in c:\\programdata\\anaconda3\\lib\\site-packages (from requests<3,>=2.27->streamlit) (2024.8.30)\n",
      "Requirement already satisfied: markdown-it-py>=2.2.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from rich<14,>=10.14.0->streamlit) (2.2.0)\n",
      "Requirement already satisfied: pygments<3.0.0,>=2.13.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from rich<14,>=10.14.0->streamlit) (2.15.1)\n",
      "Requirement already satisfied: smmap<5,>=3.0.1 in c:\\programdata\\anaconda3\\lib\\site-packages (from gitdb<5,>=4.0.1->gitpython!=3.1.19,<4,>=3.0.7->streamlit) (4.0.0)\n",
      "Requirement already satisfied: MarkupSafe>=2.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from jinja2->altair<6,>=4.0->streamlit) (2.1.3)\n",
      "Requirement already satisfied: attrs>=22.2.0 in c:\\programdata\\anaconda3\\lib\\site-packages (from jsonschema>=3.0->altair<6,>=4.0->streamlit) (23.1.0)\n",
      "Requirement already satisfied: jsonschema-specifications>=2023.03.6 in c:\\programdata\\anaconda3\\lib\\site-packages (from jsonschema>=3.0->altair<6,>=4.0->streamlit) (2023.7.1)\n",
      "Requirement already satisfied: referencing>=0.28.4 in c:\\programdata\\anaconda3\\lib\\site-packages (from jsonschema>=3.0->altair<6,>=4.0->streamlit) (0.30.2)\n",
      "Requirement already satisfied: rpds-py>=0.7.1 in c:\\programdata\\anaconda3\\lib\\site-packages (from jsonschema>=3.0->altair<6,>=4.0->streamlit) (0.10.6)\n",
      "Requirement already satisfied: mdurl~=0.1 in c:\\programdata\\anaconda3\\lib\\site-packages (from markdown-it-py>=2.2.0->rich<14,>=10.14.0->streamlit) (0.1.0)\n",
      "Requirement already satisfied: six>=1.5 in c:\\programdata\\anaconda3\\lib\\site-packages (from python-dateutil>=2.8.2->pandas<3,>=1.3.0->streamlit) (1.16.0)\n"
     ]
    }
   ],
   "source": [
    "pip install streamlit"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 108,
   "id": "ad7ecf9a-a56c-48f6-884d-2b2ebe354cc5",
   "metadata": {},
   "outputs": [
    {
     "ename": "NameError",
     "evalue": "name 'app' is not defined",
     "output_type": "error",
     "traceback": [
      "\u001b[1;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[1;31mNameError\u001b[0m                                 Traceback (most recent call last)",
      "Cell \u001b[1;32mIn[108], line 1\u001b[0m\n\u001b[1;32m----> 1\u001b[0m app\u001b[38;5;241m.\u001b[39mpy\n",
      "\u001b[1;31mNameError\u001b[0m: name 'app' is not defined"
     ]
    }
   ],
   "source": [
    "app.py"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 110,
   "id": "36bd19cb-42c8-43a8-8193-a4c96bae2422",
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "2026-03-12 13:00:45.696 \n",
      "  \u001b[33m\u001b[1mWarning:\u001b[0m to view this Streamlit app on a browser, run it with the following\n",
      "  command:\n",
      "\n",
      "    streamlit run C:\\ProgramData\\anaconda3\\Lib\\site-packages\\ipykernel_launcher.py [ARGUMENTS]\n",
      "2026-03-12 13:00:45.702 Session state does not function when running a script without `streamlit run`\n"
     ]
    }
   ],
   "source": [
    "import streamlit as st\n",
    "import joblib\n",
    "import pandas as pd\n",
    "\n",
    "st.title(\"🧬 Huntington Disease Risk Predictor\")\n",
    "\n",
    "repeat = st.number_input(\n",
    "    \"Enter CAG Repeat Count\",\n",
    "    min_value=0,\n",
    "    max_value=100,\n",
    "    value=20\n",
    ")\n",
    "\n",
    "model = joblib.load(\"huntington_ml_model.pkl\")\n",
    "\n",
    "if st.button(\"Predict\"):\n",
    "    \n",
    "    result = model.predict(\n",
    "        pd.DataFrame({\"CAG_Repeats\":[repeat]})\n",
    "    )[0]\n",
    "    \n",
    "    if result == 1:\n",
    "        st.error(\"⚠️ High Risk of Huntington Disease\")\n",
    "    else:\n",
    "        st.success(\"✅ Normal Range\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "0078a508-ce9f-450a-8825-1741eca5b569",
   "metadata": {},
   "outputs": [],
   "source": [
    "\n"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python [conda env:anaconda3]",
   "language": "python",
   "name": "conda-env-anaconda3-py"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.12.7"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
