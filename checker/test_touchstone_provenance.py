"""Constructed transport controls; these values are not Excel observations."""
import copy
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from racer import fingerprint as f


class TouchstoneProvenance(unittest.TestCase):
    def fixture(self,result_index=None):
        client,_=f.load_service(f.ROOT.parent/'DnaCalc')
        from touchstone import protocol
        rows=[[client.number(1.)]];rid='1'*32;env_id='0123456789ab'
        env=dict(engine='excel',version='16.0',build='20527',file_version='16.0.20527.20150',
                 os='constructed',cpu_id='constructed',mt_calc='off',iteration=False,
                 automation_security=1,precision_as_displayed=False,date1904=False,
                 compatibility_version=2,use_system_separators=True,decimal_separator='.',
                 thousands_separator=',',calculation=-4135)
        reply=dict(version=1,id=rid,ok=True,inputs=rows,rows=[client.number(2.)],
                   automation_security_before=1,automation_security_after=1,
                   excel_version=env['version'],excel_build=env['build'],excel_file_version=env['file_version'],
                   machine='constructed',timestamp='constructed',request_sha256=None,env_id=env_id,
                   calculation_settings=dict(PrecisionAsDisplayed=False,Date1904=False,CompatibilityVersion=2,
                                             UseSystemSeparators=True,DecimalSeparator='.',ThousandsSeparator=',',calculation=-4135))
        groups=protocol.rows_to_groups('ABS',rows,result_index)
        pending=SimpleNamespace(rid=rid,env_ids={env_id},envs={env_id:dict(env=env)},
             frames_sent=protocol.submit_frames(rid,groups,dict(echo=True,order='any',calls=1)))
        return rows,reply,pending

    def test_frame_and_environment_evidence_accepts_unchanged_reply(self):
        rows,reply,pending=self.fixture();original=copy.deepcopy(reply)
        audit=f.touchstone_audit(pending,reply)
        f.accept(reply,'ABS',rows,audit=audit)
        self.assertEqual(reply,original)

    def test_priority_is_transport_metadata_but_other_wire_fields_stay_bound(self):
        rows,reply,pending=self.fixture()
        for priority in ['interactive','bulk']:
            changed=copy.deepcopy(pending)
            changed.frames_sent[0]['opts']['priority']=priority
            f.accept(reply,'ABS',rows,audit=f.touchstone_audit(changed,reply))
        for field,value in [('priority','unknown'),('order','serial'),('calls',2)]:
            changed=copy.deepcopy(pending);changed.frames_sent[0]['opts'][field]=value
            with self.assertRaises(ValueError):f.accept(reply,'ABS',rows,audit=f.touchstone_audit(changed,reply))

    def test_environment_id_alone_never_replaces_legacy_receipt(self):
        rows,reply,_=self.fixture()
        with self.assertRaises(ValueError):f.accept(reply,'ABS',rows)
        reply['session_receipt_sha256']=None
        with self.assertRaises(ValueError):f.accept(reply,'ABS',rows)

    def test_reply_and_audit_mutations_fail(self):
        rows,reply,pending=self.fixture();audit=f.touchstone_audit(pending,reply)
        for name in ['rows','inputs','automation_security_after','excel_build','env_id']:
            changed=copy.deepcopy(reply);changed[name]=None
            with self.assertRaises(ValueError):f.accept(changed,'ABS',rows,audit=audit)
        changed=copy.deepcopy(audit);changed['submit_frames'][0]['id']='2'*32
        with self.assertRaises(ValueError):f.accept(reply,'ABS',rows,audit=changed)

    def test_valid_hashes_do_not_excuse_wrong_wire_function(self):
        rows,reply,pending=self.fixture();audit=f.touchstone_audit(pending,reply)
        with self.assertRaises(ValueError):f.accept(reply,'EXP',rows,audit=audit)

    def test_environment_build_and_settings_must_match(self):
        rows,reply,pending=self.fixture()
        for field,value in [('build','1'),('automation_security',3),('date1904',True),('precision_as_displayed',True)]:
            changed=copy.deepcopy(pending);changed.envs['0123456789ab']['env'][field]=value
            audit=f.touchstone_audit(changed,reply)
            with self.assertRaises(ValueError):f.accept(reply,'ABS',rows,audit=audit)

    def test_v2_cache_keeps_raw_reply_and_audit_apart(self):
        rows,reply,pending=self.fixture();audit=f.touchstone_audit(pending,reply)
        client=SimpleNamespace(last_audit=audit,ask=lambda *args,**kwargs:copy.deepcopy(reply))
        with tempfile.TemporaryDirectory()as directory:
            live=f.CachedService(client,directory)
            self.assertEqual(live.ask('ABS',rows),reply)
            offline=f.CachedService(None,directory,True)
            self.assertEqual(offline.ask('ABS',rows),reply)
            path=next(Path(directory).glob('*.json.gz'))
            import gzip
            with gzip.open(path,'rt')as stream:stored=json.load(stream)
            self.assertEqual(stored['reply'],reply)
            self.assertEqual(stored['audit'],audit)

    def test_mixin_requests_input_echoes(self):
        rows,reply,pending=self.fixture();received={}
        pending.result=lambda timeout:copy.deepcopy(reply)
        def submit(function,arguments,**options):
            received.update(function=function,rows=arguments,options=options)
            return pending
        connection=SimpleNamespace(submit_rows=submit)
        class Base:
            timeout=7
            def _socket(self):return connection
        class Client(f.AuditedClientMixin,Base):pass
        client=Client();actual=client.ask('ABS',rows)
        self.assertTrue(received['options']['echo'])
        self.assertEqual(actual,reply)
        f.accept(actual,'ABS',rows,audit=client.last_audit)

    def test_existing_cached_service_callers_get_audited_socket_replies(self):
        rows,reply,pending=self.fixture();options={}
        pending.result=lambda timeout:copy.deepcopy(reply)
        def submit(function,arguments,**kwargs):options.update(kwargs);return pending
        connection=SimpleNamespace(submit_rows=submit)
        class PlainClient:
            timeout=7
            def _socket(self):return connection
            def ask(self,*args,**kwargs):raise AssertionError('socket call bypassed its audit')
        with tempfile.TemporaryDirectory()as directory:
            service=f.CachedService(PlainClient(),directory)
            self.assertEqual(service.ask('ABS',rows),reply)
            self.assertTrue(options['echo'])
            self.assertEqual(f.CachedService(None,directory,True).ask('ABS',rows),reply)


if __name__=='__main__':unittest.main()
