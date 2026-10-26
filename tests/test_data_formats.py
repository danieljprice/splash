"""CLI data-format tests: splash --labels / --header and temp ascii conversion."""
from pathlib import Path
import json
import os
import re
import shutil
import subprocess

import pytest

# Formats with no small dump in splash-testdata yet (too large, hangs, or extra libs):
# dragon, mbate, jjm, spyros, sro/magma, mhutch, ucla, foulkes, vanaverbeke,
# gadget_jsb, egaburov, cactus_hdf5, phantom_hdf5, shamrock, h5part, pbob,
# starsmasher, oilonwater.

HEADER_RE = re.compile(r"^\s*(\S+)\s*=\s*(\S+)")
ABORT_CODES = {132, 134, 136, 138, 139}
LABELS_TIMEOUT = 8
CONVERT_TIMEOUT = 30


def load_format_cases(testdata_root: Path) -> dict:
    path = testdata_root / "expected_formats.json"
    assert path.is_file(), f"missing {path}"
    with path.open() as fh:
        data = json.load(fh)
    assert data.get("schema") == 1
    assert data.get("cases")
    return data


def splash_exe() -> str:
    return os.environ.get("SPLASH_BIN", "splash")


def run_splash(args, cwd=None, env_extra=None, timeout=LABELS_TIMEOUT):
    env = os.environ.copy()
    if env_extra:
        env.update(env_extra)
    proc = subprocess.run(
        [splash_exe(), *args],
        cwd=cwd,
        env=env,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        timeout=timeout,
    )
    proc.stdout = (proc.stdout or b"").decode("utf-8", errors="replace")
    proc.stderr = (proc.stderr or b"").decode("utf-8", errors="replace")
    return proc


def aborted(proc) -> bool:
    if proc.returncode in ABORT_CODES:
        return True
    blob = (proc.stdout or "") + "\n" + (proc.stderr or "")
    return "Fortran runtime error" in blob or "Backtrace for this error" in blob


def whole_lines(text: str) -> set:
    return {line.strip() for line in text.splitlines() if line.strip()}


def parse_header(text: str) -> dict:
    out = {}
    for line in text.splitlines():
        match = HEADER_RE.match(line)
        if match:
            out[match.group(1)] = match.group(2)
    return out


def _formats_text():
    proc = run_splash(["--formats"], timeout=5)
    return (proc.stdout or "") + (proc.stderr or "")


def _build_features():
    text = _formats_text()
    return {
        "hdf5": "This build does not support HDF5" not in text,
        "fits": "This build does not support FITS" not in text,
    }


def _case_by_id(testdata_root: Path, case_id: str) -> dict:
    data = load_format_cases(testdata_root)
    for case in data["cases"]:
        if case["id"] == case_id:
            return case
    raise KeyError(case_id)


def _maybe_skip(case, features):
    for req in case.get("requires") or []:
        if req == "hdf5" and not features["hdf5"]:
            pytest.skip("splash was not compiled with HDF5")
        if req == "fits" and not features["fits"]:
            pytest.skip("splash was not compiled with FITS")


def _splash_args(case, dump: Path, kind: str):
    args = []
    if kind == "labels":
        args.append("--labels")
    elif kind == "header":
        args.append("--header")
    fmt = case.get("format")
    if fmt:
        args.extend(["-f", fmt])
    args.extend(case.get("extra_args") or [])
    args.append(str(dump))
    return args


class TestDataFormats:

    def test_expectations_file(self, splash_testdata):
        data = load_format_cases(splash_testdata)
        ids = [case["id"] for case in data["cases"]]
        assert len(ids) == len(set(ids)), "duplicate case ids in expected_formats.json"
        for case in data["cases"]:
            path = splash_testdata / case["file"]
            assert path.is_file(), f"missing dump {path}"
            for extra in case.get("companion_files") or []:
                assert (splash_testdata / extra).is_file(), f"missing companion {extra}"

    @pytest.fixture(scope="class")
    def features(self):
        return _build_features()

    def test_labels(self, splash_testdata, features, case_id):
        case = _case_by_id(splash_testdata, case_id)
        _maybe_skip(case, features)
        dump = splash_testdata / case["file"]
        proc = run_splash(
            _splash_args(case, dump, "labels"),
            env_extra=case.get("env"),
        )
        assert not aborted(proc), (
            f"{case_id}: splash aborted\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
        if case.get("expect_success", True):
            assert proc.returncode == 0, (
                f"{case_id}: splash --labels failed ({proc.returncode})\n"
                f"{proc.stdout}\n{proc.stderr}"
            )
            lines = whole_lines(proc.stdout)
            for label in case.get("labels") or []:
                assert label in lines, (
                    f"{case_id}: expected label {label!r} as a whole line\n"
                    f"stdout:\n{proc.stdout}"
                )
            for snippet in case.get("stdout_contains") or []:
                blob = (proc.stdout or "") + (proc.stderr or "")
                assert snippet in blob, f"{case_id}: expected {snippet!r} in output"
        else:
            blob = (proc.stdout or "") + (proc.stderr or "")
            assert "ERROR" in blob or proc.returncode != 0

    def test_header(self, splash_testdata, features, case_id):
        case = _case_by_id(splash_testdata, case_id)
        _maybe_skip(case, features)
        expected = case.get("header") or {}
        if not expected:
            pytest.skip("no header checks for this case")
        dump = splash_testdata / case["file"]
        proc = run_splash(
            _splash_args(case, dump, "header"),
            env_extra=case.get("env"),
        )
        assert not aborted(proc), f"{case_id}: splash --header aborted"
        assert proc.returncode == 0, (
            f"{case_id}: splash --header failed ({proc.returncode})\n"
            f"{proc.stdout}\n{proc.stderr}"
        )
        got = parse_header(proc.stdout)
        for tag, want in expected.items():
            assert tag in got, f"{case_id}: missing header tag {tag}\n{proc.stdout}"
            gotv = float(got[tag])
            wantv = float(want)
            if wantv == 0:
                assert gotv == 0, f"{case_id}: {tag}={gotv}, expected 0"
            else:
                rel = abs(gotv - wantv) / abs(wantv)
                assert rel <= 1.0e-4, f"{case_id}: {tag}={gotv}, expected {wantv}"

    def test_auto_detect(self, splash_testdata, features, case_id):
        case = _case_by_id(splash_testdata, case_id)
        _maybe_skip(case, features)
        if not case.get("auto"):
            pytest.skip("auto-detect not expected for this case")
        if not case.get("expect_success", True):
            pytest.skip("robustness case")
        dump = splash_testdata / case["file"]
        args = ["--labels", *(case.get("extra_args") or []), str(dump)]
        proc = run_splash(args, env_extra=case.get("env"))
        assert not aborted(proc)
        assert proc.returncode == 0, (
            f"{case_id}: auto-detect --labels failed ({proc.returncode})\n"
            f"{proc.stdout}\n{proc.stderr}"
        )
        lines = whole_lines(proc.stdout)
        for label in case.get("labels") or []:
            assert label in lines, f"{case_id}: auto-detect missing label {label!r}"

    def test_convert_ascii(self, splash_testdata, features, tmp_path, case_id):
        case = _case_by_id(splash_testdata, case_id)
        _maybe_skip(case, features)
        if not case.get("convert_ascii"):
            pytest.skip("ascii conversion not requested")
        src = splash_testdata / case["file"]
        dump = tmp_path / src.name
        shutil.copy2(src, dump)
        for extra in case.get("companion_files") or []:
            shutil.copy2(splash_testdata / extra, tmp_path / Path(extra).name)
        args = ["to", "ascii"]
        if case.get("format"):
            args.extend(["-f", case["format"]])
        args.extend(case.get("extra_args") or [])
        args.append(src.name)
        proc = run_splash(args, cwd=tmp_path, env_extra=case.get("env"),
                          timeout=CONVERT_TIMEOUT)
        assert not aborted(proc), f"{case_id}: splash to ascii aborted"
        assert proc.returncode == 0, (
            f"{case_id}: splash to ascii failed ({proc.returncode})\n"
            f"{proc.stdout}\n{proc.stderr}"
        )
        ascii_files = [p for p in tmp_path.iterdir() if p.name.endswith(".ascii")]
        assert ascii_files, (
            f"{case_id}: no .ascii output in {tmp_path}\n"
            f"files: {list(tmp_path.iterdir())}\n"
            f"{proc.stdout}\n{proc.stderr}"
        )
        assert any(path.stat().st_size > 0 for path in ascii_files), (
            f"{case_id}: ascii output is empty"
        )


def pytest_generate_tests(metafunc):
    if "case_id" not in metafunc.fixturenames:
        return
    env = os.environ.get("SPLASH_TESTDATA")
    if not env:
        metafunc.parametrize("case_id", ["skip-no-testdata"])
        return
    data = load_format_cases(Path(env))
    metafunc.parametrize("case_id", [case["id"] for case in data["cases"]])
