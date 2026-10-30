from pathlib import Path
import pytest
from wmh_gmv import launcher


@pytest.mark.parametrize("action,module", [("start","wmh_gmv.lifecycle"),("stop","wmh_gmv.lifecycle"),("monitor","wmh_gmv.monitor")])
def test_configured_interpreter_is_authoritative(tmp_path,monkeypatch,action,module):
    expected=str(tmp_path/'custom environment/bin/python')
    config={"paths":{"python":expected},"project":{"project_root":str(tmp_path)},"_config_path":str(tmp_path/'profile.yaml')}
    monkeypatch.setattr(launcher,"load_config",lambda path:config)
    captured={}
    def execute(executable,args,env):
        captured.update(executable=executable,args=args,env=env)
        raise RuntimeError('intercept exec')
    monkeypatch.setattr(launcher.os,"execve",execute)
    with pytest.raises(RuntimeError,match='intercept exec'):
        launcher.main([action,'--config','profile.yaml',*(['--once'] if action=='monitor' else [])])
    assert captured['executable']==expected
    assert captured['args'][:3]==[expected,'-m',module]
    assert str(tmp_path/'src') in captured['env']['PYTHONPATH']
    assert captured['args'][-1]=='--once' if action=='monitor' else captured['args'][-1]==config['_config_path']
