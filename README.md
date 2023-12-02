# CA_final_project
Spec: 
    ◆All inputs are synchronized with the negative edge clock.
    ◆All outputs should be synchronized at clock rising edge.
    ◆You should reset all your outputs when i_rst_n is low.
Todo:
    *Submodule
    1. ALU
        (i: 2 data o: 1 data)
        
    2. Control unit
        (i: opcode o: control signals)
        input opcode to determine which control signal need to set
    
    *always block
    1. Intruction fetch
        Set high o_IMEM_cen to load instruction
        pass PC address
        return i_IMEM_data(instruction)

    2. nextPC <= PC+4

Q:
    load store有沒有要signed extention?
    Mux怎麼放?
    branch或ALU control是弄一個always block嗎?
