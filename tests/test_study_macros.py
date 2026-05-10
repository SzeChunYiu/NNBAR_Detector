from pathlib import Path


def test_pi0_foil_energy_scan_macro_runs_fixed_energy_ladder_from_foil() -> None:
    macro = Path("macro/studies/pi0_foil_energy_scan.mac")

    text = macro.read_text(encoding="utf-8")

    assert "/calibration/mode true" not in text
    assert "/particle_generator/set_folder_name studies/pi0_foil_energy_scan" in text
    assert "/calibration/signal_particle pi0" in text
    for run, energy in enumerate((100, 200, 300, 400, 500, 600)):
        assert f"/calibration/signal_energy_min {energy} MeV" in text
        assert f"/calibration/signal_energy_max {energy} MeV" in text
        assert f"/particle_generator/set_run_number {run}" in text
    assert text.count("/particle_generator/set_event_number 0") == 6
    assert text.count("/run/beamOn 200") == 6


def test_charged_pion_proton_stress_macro_uses_signal_mode_from_foil() -> None:
    macro = Path("macro/studies/charged_pion_proton_foil_stress.mac")

    text = macro.read_text(encoding="utf-8")

    assert "/calibration/mode true" not in text
    assert "/particle_generator/set_folder_name studies/charged_pion_proton_foil_stress" in text
    assert "/calibration/signal_particle pi+" in text
    assert "/calibration/signal_particle pi-" in text
    assert "/calibration/signal_particle proton" in text
    assert "/calibration/signal_energy_min 50 MeV" in text
    assert "/calibration/signal_energy_max 600 MeV" in text
    assert "/calibration/signal_energy_max 500 MeV" in text
    assert "/particle_generator/set_run_number 0" in text
    assert "/particle_generator/set_run_number 1" in text
    assert "/particle_generator/set_run_number 2" in text
    assert text.count("/run/beamOn 200") == 3


def test_multiprimary_pion_proton_stress_macro_emits_shared_event_topology() -> None:
    macro = Path("macro/studies/multiprimary_pion_proton_foil_stress.mac")

    text = macro.read_text(encoding="utf-8")

    assert "/calibration/mode true" not in text
    assert "/particle_generator/set_folder_name studies/multiprimary_pion_proton_foil_stress" in text
    assert "/calibration/signal_particles pi+,pi-,proton" in text
    assert "/calibration/signal_particle " not in text
    assert "/calibration/signal_energy_min 50 MeV" in text
    assert "/calibration/signal_energy_max 500 MeV" in text
    assert "/particle_generator/set_run_number 0" in text
    assert text.count("/run/beamOn 200") == 1
