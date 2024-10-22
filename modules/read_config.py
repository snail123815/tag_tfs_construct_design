from pathlib import Path
import json
import logging

logger = logging.getLogger(__name__)

def read_config(config_file:Path) -> dict:
    """
    Read a config file and return a dictionary of the config.
    Checks essential keys are present and correct.
    """
    config = {}
    with config_file.open('r') as f:
        json_config = json.load(f)
        for entry in json_config:
            config[entry] = json_config[entry]
    logger.info(f"Read config file: {config_file}")

    assert "hmmsearch" in config, "hmmsearch path not found in config"
    assert Path(config["hmmsearch"]).is_file(), "hmmsearch path not found"
    assert "gather_threshold_e" in config, "gather_threshold_e not found in config"
    assert "domain_coverage" in config, "domain_coverage not found in config"

    config["hmmsearch"] = Path(config["hmmsearch"])
    config["gather_threshold_e"] = float(config["gather_threshold_e"])
    config["domain_coverage"] = float(config["domain_coverage"])

    return config