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
    assert config["hmmsearch"].is_file(), "hmmsearch path not found"

    return config