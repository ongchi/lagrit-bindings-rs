import pytest
import doctest


@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    output_path = request.path.parent / "pytest_output"
    output_path.mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(output_path)


def test_pylagrit():
    doctest.testfile("../pylagrit/core.py", verbose=False, raise_on_error=True)
