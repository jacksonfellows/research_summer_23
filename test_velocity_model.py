import velocity_model
import profile_info

# Make sure we can load test_01.vm.


def test_load_test_01():
    vm = velocity_model.VMTOMO_VM.load("test_data/test_01.vm")
    assert vm.x1 == 0.0
    assert vm.x2 == 250.0
    assert vm.z1 == -2.0
    assert vm.z2 == 60.0
    assert vm.nx == 750
    assert vm.nz == 300
    assert vm.nr == 3


def test_plot_test_01():
    vm = velocity_model.VMTOMO_VM.load("test_data/test_01.vm")
    vm.plot(show=False)  # See if it runs without errors.


def test_build_profile():
    p1 = velocity_model.build_vm_1()
    p1.plot(show=False)
