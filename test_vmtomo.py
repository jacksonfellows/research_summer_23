import vmtomo

# Make sure we can load test_01.vm.


def test_load_test_01():
    vm = vmtomo.VMTOMO_VM("test_01.vm")
    assert vm.x1 == 0.0
    assert vm.x2 == 250.0
    assert vm.z1 == -2.0
    assert vm.z2 == 60.0
    assert vm.nx == 750
    assert vm.nz == 300
    assert vm.nr == 3
