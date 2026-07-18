from plinder.data import databases


def test_make_db_limits_external_tool_threads(tmp_path, monkeypatch):
    commands: list[list[str]] = []
    monkeypatch.setattr(databases, "run", commands.append)

    databases.make_db(
        input_dir=tmp_path / "input",
        output_dir=tmp_path / "output",
        db="foldseek",
        threads=3,
    )

    assert len(commands) == 2
    for command in commands:
        index = command.index("--threads")
        assert command[index + 1] == "3"
