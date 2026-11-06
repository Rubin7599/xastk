import json

import numpy as np

from xastk import __version__
from xastk.cli import main


def test_nidp(capsys):
    assert main(["nidp", "3", "17", "1.4", "3.0"]) == 0
    assert capsys.readouterr().out.strip() == "15.260"


def test_hamilton(capsys):
    main(["hamilton", "--chi2-a", "112.7", "--nvar-a", "4", "--chi2-b", "68.3", "--nvar-b", "5", "--nidp", "15.26"])
    assert "significant at 95%" in capsys.readouterr().out


def test_table_and_check(synthetic_result, tmp_path, capsys):
    rpath = tmp_path / "r.json"
    rpath.write_text(json.dumps(synthetic_result))
    main(["table", str(rpath), "--to", "csv", "latex", "--out", str(tmp_path / "si")])
    assert (tmp_path / "si.csv").exists() and (tmp_path / "si.tex").exists()
    assert main(["check", str(rpath)]) == 0
    assert main(["check", str(rpath), "--strict"]) == 1
    assert "correlation" in capsys.readouterr().out


def test_wt_command(tmp_path, capsys):
    import matplotlib

    matplotlib.use("Agg")
    k = np.arange(0, 15, 0.05)
    np.savetxt(tmp_path / "chi.dat", np.column_stack([k, np.sin(2 * k * 2.5) / np.clip(k, 0.5, None)]))
    main(["wt", str(tmp_path / "chi.dat"), "--out", str(tmp_path / "wt.png")])
    assert "lobe:" in capsys.readouterr().out
    assert (tmp_path / "wt.png").exists()


def test_version():
    assert __version__.startswith("0.1")
