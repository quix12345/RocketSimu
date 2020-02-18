class nozzle:
    def __init__(self, nozzle_type="straight", effiency=0.85, aeat=1):
        self.nozzle_type = nozzle_type
        self.effiency = effiency
        if self.nozzle_type == "customize_laval":
            self.aeat = aeat
        elif self.nozzle_type == "auto_laval":
            self.aeat = -1
        elif self.nozzle_type == 'straight':
            self.aeat = 1
        else:
            raise Exception("Unknown nozzle type!",nozzle_type)

class straight_nozzle(nozzle):
    def __init__(self, effiency=0.85):
        super(straight_nozzle, self).__init__(effiency=effiency)

class laval_nozzle(nozzle):
    def __init__(self, effiency=0.85, IsCustomized=False,aeat=1):
        if IsCustomized:
            nozzle_type="customize_laval"
            super(laval_nozzle, self).__init__(nozzle_type=nozzle_type, effiency=effiency,aeat=aeat)
        else:
            nozzle_type= "auto_laval"
            super(laval_nozzle, self).__init__(nozzle_type=nozzle_type, effiency=effiency)
